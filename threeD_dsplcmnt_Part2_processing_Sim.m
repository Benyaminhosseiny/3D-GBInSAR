% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

% * First, run "threeD_dsplcmnt_Part1_rawdata.m"
% * Or
% * import your own "cube_3dTS" (a 4D array of TS of raw cube data [r x az x el x time])

% clear;clc;close all
c = physconst('LightSpeed');

%% Input setting:
% % SAR Processing:
r_fft  = size(cube_3dTS,1);
az_fft = 64;
el_fft = 64;

% % InSAR Processing:
adi_thresh = 0.2;
amp_thresh = 2; % dB from the max

% =========================================================================

%% Known Parameters:
Nr = size(cube_3dTS,1);
Na = size(cube_3dTS,2);
Ne = size(cube_3dTS,3);
Nts= size(cube_3dTS,4); % Number of TS data
d_az = d_az;
d_el = d_el;

beam_az  = rad2deg(lambda/2/d_az);
beam_el  = rad2deg(lambda/2/d_el);
R_ax     = linspace(rr/2,rr*Nr-rr/2,r_fft)'; 
Theta_ax = -linspace(-beam_az/2,beam_az/2,az_fft);
Phi_ax   = -linspace(-beam_el/2,beam_el/2,el_fft);

% Radar antennas pos:
X_rad = d_az*[-floor(Na/2):floor(Na/2)];
Y_rad0 = 0;
Z_rad0 = d_el*[-floor(Ne/2):floor(Ne/2)];
nadir=-0;
rot_mat=[1 0 0; 0 cosd(nadir) -sind(nadir); 0 sind(nadir) cosd(nadir)];
z_shift=0;

Y_rad = Y_rad0*cosd(nadir)-Z_rad0*sind(nadir);
Z_rad = z_shift+Y_rad0*sind(nadir)+Z_rad0*cosd(nadir);

%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================

%% SAR Processing:
for ts_ii = 1:Nts
    cube_3d=cube_3dTS(:,:,:,ts_ii);
% % Windowing
%     cube_3d=windowing(cube_3d,[1,2,3]);
    
% % % SLC
    [rc, slc3d] = create_slc_3d(cube_3d,r_fft,az_fft,el_fft);
    slc3d = slc3d./max(abs(slc3d),[],'all');
    
    slc3dTS(:,:,:,ts_ii) = slc3d(:,:,:);  % Fully-Compressed (3D SLC) time-series data-cube
    rcTS   (:,:,:,ts_ii) = rc(:,:,:);     % Range-Compressed time-series data-cube
end

slc_db = 10*log10(abs(slc3d)); 
slc_db = slc_db-max(slc_db,[],'all'); % Set maximum backscatter to zero
clear cube_3dTS

%% InSAR Processing:
% % % 3D ADI analysis:
adi = ADI_3d(slc3dTS);

% % % PS detection:
ps_mask = PS_3d(slc3dTS, adi_thresh, amp_thresh);%ps_mask(ps_mask==0)=NaN;

% % % point cloud:
point_cloud = gen_pointcloud(ps_mask(:,:,:),0);
ps_idxr=point_cloud(:,1);ps_idxaz=point_cloud(:,2);ps_idxel=point_cloud(:,3);


for ii=1:length(ps_idxr)
    ps_mag(ii,1)  = slc_db(ps_idxr(ii),ps_idxaz(ii),ps_idxel(ii));
    ps_adi(ii,1)  = adi(ps_idxr(ii),ps_idxaz(ii),ps_idxel(ii));
    ps_R(ii,1)    = R_ax(ps_idxr(ii));
    ps_Theta(ii,1)= Theta_ax(ps_idxaz(ii));
    ps_Phi(ii,1)  = Phi_ax(ps_idxel(ii));
end

% Cartesian:
ps_X2=ps_R.*cosd(ps_Phi).*sind(ps_Theta);
ps_Y2=ps_R.*cosd(ps_Phi).*cosd(ps_Theta);
ps_Z2=ps_R.*sind(ps_Phi);

ps_X=ps_R.*sind(ps_Theta);
ps_Z=ps_R.*sind(ps_Phi);
ps_Y=ps_R.*sqrt( 1-(sind(ps_Theta).^2+sind(ps_Phi).^2) );

% ps_X=X_tar*ones(size(ps_Theta));
% ps_Z=Z_tar*ones(size(ps_Theta));
% ps_Y=Y_tar*ones(size(ps_Theta));


% % % LOS:
[TS_phase,TS_intf_phase,TS_cum_phase] = TSInSAR_3d(slc3dTS);
for ii=1:length(ps_idxr)
    ps_los(ii)=TS_cum_phase(ps_idxr(ii),ps_idxaz(ii),ps_idxel(ii),end)*1e3*lambda/4/pi;
    ps_losTS(:,ii)=squeeze(TS_cum_phase(ps_idxr(ii),ps_idxaz(ii),ps_idxel(ii),:))'*1e3*lambda/4/pi;
end

% % % 3D:
% % % % Detected PS' XYZ:
Tar3dLoc=[ps_X,ps_Y,ps_Z];
Tar3dLoc=(rot_mat*Tar3dLoc')';
Tar3dLoc(:,3)=z_shift+Tar3dLoc(:,3);
ps_X=Tar3dLoc(:,1);
ps_Y=Tar3dLoc(:,2);
ps_Z=Tar3dLoc(:,3);
% % % % Sensor steps' XYZ:
x_sensor =  repmat(X_rad', 1,length(Z_rad));
y_sensor =  repmat(Y_rad, length(X_rad),1);
z_sensor =  repmat(Z_rad, length(X_rad),1);

% Radar's 3D position during data aquisition for SAR imaging. 
%[3d array NxMx3 ]:
SAR3dLoc(:,:,1)=x_sensor;
SAR3dLoc(:,:,2)=y_sensor;
SAR3dLoc(:,:,3)=z_sensor;

clear RC_sig_PS_ts

% TS RC signals corresponding to the detected PS:
for ii=1:length(ps_idxr)
    rcii=rcTS(ps_idxr(ii),:,:,:);
    rcii=reshape(rcii,Na*Ne,Nts);
    RC_sig_PS_ts(ii,:,:)=rcii;
end
RC_sig_PS_ts=permute(RC_sig_PS_ts,[3,1,2]);

clutter_rmv_flag=0;
% % %
[d_hat_3D_cartesian_dif, d_hat_3D_cartesian_ts, d_hat_3D_cartesian_total] = Displacement_vec_Cartesian(Tar3dLoc, SAR3dLoc, RC_sig_PS_ts, lambda, clutter_rmv_flag);
ps_dX   = d_hat_3D_cartesian_total(:,1)*1e3;
ps_dY   = d_hat_3D_cartesian_total(:,2)*1e3;
ps_dZ   = d_hat_3D_cartesian_total(:,3)*1e3;
ps_dXTS = squeeze(d_hat_3D_cartesian_ts(:,:,1))*1e3; %t*p mm
ps_dYTS = squeeze(d_hat_3D_cartesian_ts(:,:,2))*1e3; %t*p mm
ps_dZTS = squeeze(d_hat_3D_cartesian_ts(:,:,3))*1e3; %t*p mm

%% Show point clouds:
fontsizefig = 14; fontname = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',fontsizefig);

xlimit=[-2+min(X_tar), 2+max(X_tar)];
ylimit=[-5+min(Y_tar), 5+max(Y_tar)];
zlimit=[-2+min(Z_tar), 2+max(Z_tar)];
figure('Position', [20 100 1500 900]); sgtitle("3D Scatterers (axes are in [m])")
subplot(2,3,1); scatter3_PC(ps_Y,ps_X,ps_Z,ps_mag, 'jet', 'amplitude'                  , "[dB]", xlimit, ylimit, zlimit); 
hold on; scatter3(Y_tar,X_tar,Z_tar,400,'r');
subplot(2,3,2); scatter3_PC(ps_Y,ps_X,ps_Z,ps_los, 'jet', 'LOS Cumulative displacement', "[mm]", xlimit, ylimit, zlimit); 
subplot(2,3,3); scatter3_PC(ps_Y,ps_X,ps_Z,ps_adi, 'jet', 'ADI'                        , ""    , xlimit, ylimit, zlimit); 
subplot(2,3,4); scatter3_PC(ps_Y,ps_X,ps_Z,ps_dX , 'jet', 'X Cumulative displacement'  , "[mm]", xlimit, ylimit, zlimit); 
subplot(2,3,5); scatter3_PC(ps_Y,ps_X,ps_Z,ps_dY , 'jet', 'Y Cumulative displacement'  , "[mm]", xlimit, ylimit, zlimit);
subplot(2,3,6); scatter3_PC(ps_Y,ps_X,ps_Z,ps_dZ , 'jet', 'Z Cumulative displacement'  , "[mm]", xlimit, ylimit, zlimit);

xlimit=[-5+min(ps_idxaz), 5+max(ps_idxaz)];
ylimit=[-5+min(ps_idxr),  5+max(ps_idxr)];
zlimit=[-5+min(ps_idxel), 5+max(ps_idxel)];
figure('Position', [20 100 1500 700]); sgtitle("3D Scatterers (axes are in [samples])")
subplot(2,3,1); scatter3_PC(ps_idxr,ps_idxaz,ps_idxel,ps_mag, 'jet', 'amplitude'                  , "[dB]", xlimit, ylimit, zlimit);
subplot(2,3,2); scatter3_PC(ps_idxr,ps_idxaz,ps_idxel,ps_los, 'jet', 'LOS Cumulative displacement', "[mm]", xlimit, ylimit, zlimit)
subplot(2,3,3); scatter3_PC(ps_idxr,ps_idxaz,ps_idxel,ps_adi, 'jet', 'ADI'                        , ""    , xlimit, ylimit, zlimit)
subplot(2,3,4); scatter3_PC(ps_idxr,ps_idxaz,ps_idxel,ps_dX , 'jet', 'X Cumulative displacement'  , "[mm]", xlimit, ylimit, zlimit)
subplot(2,3,5); scatter3_PC(ps_idxr,ps_idxaz,ps_idxel,ps_dY , 'jet', 'Y Cumulative displacement'  , "[mm]", xlimit, ylimit, zlimit)
subplot(2,3,6); scatter3_PC(ps_idxr,ps_idxaz,ps_idxel,ps_dZ , 'jet', 'Z Cumulative displacement'  , "[mm]", xlimit, ylimit, zlimit)

%% Show TS profile:

figure('Position', [20 100 1500 900]); sgtitle("Displacement through time (reference values are shown with red circles)")
subplot(2,2,1);plot(ps_dXTS,'--^','LineWidth',1);
hold on; scatter( [1:Nts-1], 1e3*[1:Nts-1]'*dX_tar,150,'r' ) ;xlabel('time samples');ylabel('dX [mm]')
subplot(2,2,2);plot(ps_dYTS,'--^','LineWidth',1);
hold on; scatter( [1:Nts-1], 1e3*[1:Nts-1]'*dY_tar,150,'r' ); xlabel('time samples');ylabel('dY [mm]')
subplot(2,2,3);plot(ps_dZTS,'--^','LineWidth',1);
hold on; scatter( [1:Nts-1], 1e3*[1:Nts-1]'*dZ_tar,150,'r' ); xlabel('time samples');ylabel('dZ [mm]')
subplot(2,2,4);plot(ps_losTS,'--^','LineWidth',1);            xlabel('time samples');ylabel('LOS [mm]')

% RMSE
if num_tar==1 % For more than 1 scatterers manual analysis is recommended!
    RMSEX=sqrt( mean((1e3*dX_tar*[1:Nts-1]'-mean(ps_dXTS,2)).^2) )
    RMSEY=sqrt( mean((1e3*dY_tar*[1:Nts-1]'-mean(ps_dYTS,2)).^2) )
    RMSEZ=sqrt( mean((1e3*dZ_tar*[1:Nts-1]'-mean(ps_dZTS,2)).^2) )
end

% ==============================================================================================
function scatter3_PC(ps_Y,ps_X,ps_Z,intnst,cmap,fig_title,cb_label,xlimit,ylimit,zlimit)
scatter3(ps_Y,ps_X,ps_Z,[],intnst,'filled');
axis auto
xlim(ylimit);ylim(xlimit);zlim(zlimit)
colormap(cmap);
title(fig_title);
cb=colorbar();cb.Title.String = cb_label;


end