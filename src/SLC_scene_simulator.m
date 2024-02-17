function [SNR, slc] = SLC_scene_simulator(SysSpec, TarLoc, BG_rcs, CR_rcs, CR_loc, noise_R, plot_results)
% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144
if nargin<7
    plot_results=0;
end
% SysSpec
% TarLoc
% R
% BG_rcs % m2
% CR_rcs % m2
% amp_noise % percentage(%)
% R_noise (or phase noise) % m

R = sqrt( TarLoc.x.^2 + TarLoc.y.^2 + TarLoc.z.^2 );

% cos_theta = TarLoc.y./R;
% rcs = BG_rcs*(cos_theta.^2); % m2; background's RCS

%%% Modelling as a combination of plate reflectors
% https://www.sciencedirect.com/topics/engineering/corner-reflector
zz1 = [SysSpec.x_n(1,1),  SysSpec.y_n(1), SysSpec.z_n(1,1)];
zz2 = [SysSpec.x_n(1,end), SysSpec.y_n(1),SysSpec.z_n(1,end)];
zz3 = [SysSpec.x_n(end,1), SysSpec.y_n(1), SysSpec.z_n(end,1)];
sensor_norm_vec = cross( zz1-zz2,zz1-zz3 ); % Outer product
sensor_norm_vec = sensor_norm_vec/norm(sensor_norm_vec);
sensor_norm_vec = reshape(sensor_norm_vec,[1,1,3]);

sensor_center_x = SysSpec.x_n(round(end/2),round(end/2));
sensor_center_y = SysSpec.y_n(round(end/2),round(end/2));
sensor_center_z = SysSpec.z_n(round(end/2),round(end/2));

vec_targets_x = TarLoc.x-sensor_center_x;
vec_targets_y = TarLoc.y-sensor_center_y;
vec_targets_z = TarLoc.z-sensor_center_z;

vec_targets(:,:,1) = vec_targets_x;
vec_targets(:,:,2) = vec_targets_y;
vec_targets(:,:,3) = vec_targets_z;

% % % Rectangles' normal vector
vec_targets_norm = sum(vec_targets.^2,3).^.5;
[Nx,Ny,Nz] = surfnorm(vec_targets_x,vec_targets_y,vec_targets_z);
Nt(:,:,1) = Nx; Nt(:,:,2) = Ny; Nt(:,:,3) = Nz;
% % % Targets' elevation and horizontal angles => antenna pattern => Pt
if plot_results
    figure("Position", [0 0 1000 1000]);
    subplot(3,2,1)
    surf(vec_targets_x,vec_targets_y,vec_targets_z,'FaceColor','cyan','FaceAlpha',0.8);
    axis equal; hold on
    quiver3(vec_targets_x,vec_targets_y,vec_targets_z, Nx,Ny,Nz,'color','r');
    title('imaging scene with normal vectors')
    
    subplot(3,2,2)
    idxr1=round(.4*size(vec_targets_x,1)); idxr2=round(.6*size(vec_targets_x,1));
    idxc1=round(.4*size(vec_targets_x,2)); idxc2=round(.6*size(vec_targets_x,2));
    surf(vec_targets_x(idxr1:idxr2,idxc1:idxc2),vec_targets_y(idxr1:idxr2,idxc1:idxc2),vec_targets_z(idxr1:idxr2,idxc1:idxc2),'FaceColor','cyan','FaceAlpha',0.8);
    axis equal; hold on
    quiver3(vec_targets_x(idxr1:idxr2,idxc1:idxc2),vec_targets_y(idxr1:idxr2,idxc1:idxc2),vec_targets_z(idxr1:idxr2,idxc1:idxc2), Nx(idxr1:idxr2,idxc1:idxc2),Ny(idxr1:idxr2,idxc1:idxc2),Nz(idxr1:idxr2,idxc1:idxc2),'color','r');
    title('small section of imaging scene with normal vectors')
end
tar_H_angle = atand(vec_targets_x./vec_targets_y);
tar_E_angle = asind(vec_targets_z./vec_targets_norm);
pattern_H_idx = round( (tar_H_angle/90)*length(SysSpec.pat_az)/2+(length(SysSpec.pat_az)/2) );
pattern_E_idx = round( (tar_E_angle/90)*length(SysSpec.pat_el)/2+(length(SysSpec.pat_el)/2) );
pat_az = SysSpec.pat_az-max(SysSpec.pat_az); % H-plane (or azimuth) antenna pattern (dB)
pat_el = SysSpec.pat_el-max(SysSpec.pat_el); % E-plane (or elevation) antenna pattern (dB)
pat_tar = pat_el(pattern_E_idx)+pat_az(pattern_H_idx); % Targets received power based on the antenna pattern
Pt_tar = SysSpec.Pt * 10.^( pat_tar/10 ); % transmitted power from sensor to the target's location
if plot_results
    subplot(3,2,3); imagesc(Pt_tar)
    title("transmitted power from sensor to the target's location")
end
% % % Calculate theta angle => RCS model formula
% theta = acos( sum(vec_targets.*sensor_norm_vec,3)./vec_targets_norm ); % acos(inner product) radian
theta = acos( sum(vec_targets.*Nt,3)./vec_targets_norm ); % acos(inner product) radian
if plot_results
    subplot(3,2,4); imagesc(theta)
    title("angle between a ractangle's normal vector and its location")
end

% % % RCS model
k = 2*pi/SysSpec.lambda; % Wavenumber
A = SysSpec.dx*SysSpec.dz;
rcs = BG_rcs*(( sin( k*sqrt(A).*sin(theta) )./( k*sqrt(A).*sin(theta) ) ).^2);
if ~isempty(find(theta==0))
    [zzr, zzc] = find(theta==0); % Maximum RCS
    rcs( zzr, zzc ) = BG_rcs;
end
if plot_results
    subplot(3,2,5); imagesc(rcs)
    title("RCS without CR")
end

% rcs( CR_loc(1),CR_loc(2) ) = CR_rcs; % CR
for ii = 1:size(CR_loc,1)
    rcs( CR_loc(ii,1),CR_loc(ii,2) ) = rcs( CR_loc(ii,1),CR_loc(ii,2) )+CR_rcs; % CR
end
if plot_results
    subplot(3,2,6); imagesc(rcs)
    title("RCS with CR")
end

% % % Calculate received SNR
SNR = zeros(size(rcs));
for r_i = 1:size(rcs,1)
    for c_i = 1:size(rcs,2)
        Pt_i = Pt_tar(r_i,c_i);
        %         Pt_i = SysSpec.Pt;
        SNR(r_i,c_i) = radareqsnr( SysSpec.lambda, R(r_i,c_i), Pt_i, SysSpec.tm, ...
            'RCS',rcs(r_i,c_i), 'Gain',SysSpec.Gain, 'Ts',SysSpec.Ts, 'Loss',SysSpec.Ls ); % NOTE that Pt should be based on the antenna's V and H patterns
    end
end
SNR = 10.^(SNR/10); % multiply by Noise to obtain the received Power!
slc = SNR.*exp( -1i*4*pi*R/SysSpec.lambda );

% % Adding Phase noise (range measurement!) to the signal
R_noisy = R+noise_R*2*( rand(size(slc))-0.5 );
slc = SNR.*exp( -1i*4*pi*R_noisy/SysSpec.lambda );


end