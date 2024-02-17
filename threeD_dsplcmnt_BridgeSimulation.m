% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

warning('off','all')
clear; clc; close all
addpath('./src')

c = physconst('LightSpeed');

fontsizefig = 14; fontname = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',fontsizefig);

aux_figs_flag = 1;
export_flag = 1;
export_directory = '\3D-GBInSAR\results\';
export_directory = [export_directory,datestr(now,'yyyymmdd_HHMM'),'\'];
if export_flag
    mkdir(export_directory)
end
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
%% SAR SYSTEM PROPERTIES:

%% 1- System specifics
SysSpec            = struct;
SysSpec.fc         = 77e9;       fc = SysSpec.fc;
SysSpec.tm         = 60e-6;      tm = SysSpec.tm;
SysSpec.lambda     = c/fc;       lambda = SysSpec.lambda;
SysSpec.Pt         = 0.0015;     Pt = SysSpec.Pt;
SysSpec.Gain_Tx    = 30;         Gain_Tx = SysSpec.Gain_Tx;
SysSpec.Gain_Rx    = 30;         Gain_Rx = SysSpec.Gain_Rx;
SysSpec.Gain       = [Gain_Tx, Gain_Rx];
SysSpec.peak_power = 1.5e-3;
Fn                 = 15; % noise figure (dB)
SysSpec.Fn         = 10^(Fn/10); Fn=SysSpec.Fn;
SysSpec.Ls         = 0;          Ls = SysSpec.Ls; % system loss (dB)
SysSpec.Ts         = 290*Fn;     Ts = SysSpec.Ts; % System temperature (K)
SysSpec.Bn         = 1/tm;       Bn = SysSpec.Bn; % Noise bandwidth
k                  = physconst('Boltzmann');
SysSpec.N          = k*Ts*Bn;    N = SysSpec.N; % System noise
SysSpec.fs         = 1*4.267e6;  fs= SysSpec.fs;
SysSpec.bw         = 2e9;        bw=SysSpec.bw;

%% 2-Noises
% Adding noise to the sensor's movement (Mechanical rail!)
noise_rail = 1e-5;% 1e-4; % m
% Adding Phase noise (range measurement!) to the signal
noise_R    = 1e-5; % 1e-5;% m

%% 3- Sensor movement locations
L_x =  0.2;%0.6; % GBSAR aperture length in azimuth
L_z =  0.2;%0.3; % GBSAR aperture length in elevation
d_CR = lambda/2; % Cross range spacing

% Sensor's coordinates at each step:
x_sensor =  -L_x/2+lambda/4 : d_CR : L_x/2; % 0.9 m
z_sensor =  -L_z/2+lambda/4 : d_CR : L_z/2; % 0.5 m
y_sensor =  0;

% Adding mechanical noise:
x_sensor =  x_sensor + noise_rail*2*( rand(size(x_sensor))-0.5 );
z_sensor =  z_sensor + noise_rail*2*( rand(size(z_sensor))-0.5 );
x_sensor =  repmat(x_sensor', 1,length(z_sensor));
z_sensor =  repmat(z_sensor, size(x_sensor,1),1);

%
SysSpec.x_n  =  x_sensor;
SysSpec.y_n  =  y_sensor;
SysSpec.z_n  =  z_sensor;

%% 4- Antenna Elevation and Horizontal patterns
% pat_az = xlsread('D:\Surveying\Papers\2022 - mmWave TI MIMO for SHM-JAG\H-plane pattern.xlsx');
% pat_az = interp1( pat_az(:,1),pat_az(:,2),-90:0.1:90,'spline' );
% phase
% pat_el = xlsread('D:\Surveying\Papers\2022 - mmWave TI MIMO for SHM-JAG\E-plane pattern.xlsx');
% pat_el = interp1( pat_el(:,1),pat_el(:,2),-90:0.1:90,'spline' );

varray_az = phased.ULA( 2,d_CR );
pat_az = pattern(varray_az,fc,-90:0.1:90,0,'Type','powerdb');
varray_el = phased.ULA( 4,d_CR );
pat_el = pattern(varray_el,fc,-90:0.1:90,0,'Type','powerdb');

SysSpec.pat_az = pat_az;
SysSpec.pat_el = pat_el;
if aux_figs_flag
    figure;
    plot(pat_az,'b'); hold on; plot(pat_el,'r'); title('Antenna patterns'); legend('H-plane','E-plane')
end
%% 5- Resolutions
res_r    =  .15; SysSpec.dr = res_r;% range res (affects the Quantized range image!)
res_x    =  2*0.2; SysSpec.dx = res_x;% azimuth spacing
res_z    =  0.2; SysSpec.dz = res_z;% elevation spacing

% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
%% IMAGING SCENE (TARGETS) PROPERTIES:
%% 1- Targets' location [Before Displacement]
%%% Consider a bridge with 30 m width {x = -15:15} and 6 m height {z = -3:3}.
%%% Consider range and cross-range resolutions are 0.5 m: dy=dx=dz=0.5.
%%% Consider radar is located at the distance of 20 m from the bridge.
%%% Hence bridge's y is equal to 20 (y = 20).
tar_BG_X =  25;
tar_BG_Z =  15;
tar_BG_Y =  20.21;

% Target's X coordinates
x1       =  -tar_BG_X/2:res_x:tar_BG_X/2-res_x; % 15 m in x direction

% Target's Z coordinates
z1       =  tar_BG_Z/2-res_z:-res_z:-tar_BG_Z/2;  % 6 m in z direction

[x1,z1]  =  meshgrid(x1,z1);
% Target's Y coordinates
y1       =  tar_BG_Y*ones(size(x1)); % Targets are located on a surface with e.g. y=40 m

% Adding noise to coordinates (TAKHALKHOL!!)
x1       = ( res_x/2*rand(size(x1)) )+x1;
y1       = ( res_r/2*rand(size(y1)) )+y1;
z1       = ( res_z/2*rand(size(z1)) )+z1;

TarLoc   = struct;
TarLoc.x = x1; TarLoc.y = y1; TarLoc.z = z1;


range_migration = abs(sqrt((max(x_sensor(:))-x1).^2+(max(y_sensor(:))-y1).^2+(max(z_sensor(:))-z1).^2)-sqrt((min(x_sensor(:))-x1).^2+(min(y_sensor(:))-y1).^2+(min(z_sensor(:))-z1).^2));
max(max(range_migration))

% Discrete locations based on the sensor sampling:
X_tar_hat = res_x*round(x1/res_x);
Y_tar_hat = res_r*round(y1/res_r);
Z_tar_hat = res_z*round(z1/res_z);% % Target Range with sensor's locations on the synthetic aperture

TarLoc_hat   = struct;
TarLoc_hat.x = x1;
TarLoc_hat.y = y1;
TarLoc_hat.z = z1;

% Bridge mask:
bridge_mask = zeros( size(x1) );
bridge_mask( 3*size(x1,1)/8:4*size(x1,1)/8, 1*size(x1,2)/11:10*size(x1,2)/11 ) = 1;
bridge_mask( size(x1,1)/2:end, 3*size(x1,2)/11:4*size(x1,2)/11 ) = 1;
bridge_mask( size(x1,1)/2:end, 5*size(x1,2)/11:6*size(x1,2)/11 ) = 1;
bridge_mask( size(x1,1)/2:end, 7*size(x1,2)/11:8*size(x1,2)/11 ) = 1;

bridge_mask_slc = bridge_mask;
bridge_mask_slc(bridge_mask_slc==0)=10e-30;
% figure;imagesc(bridge_mask)

%% 2- Background Targets' range before displacement
R1_map = sqrt( (x1-0).^2 + (y1-0).^2 + (z1-0).^2 ); % Before displacement (1st epoch)
% R1_map = bridge_mask .* R1_map;
% cos_theta = y./R1_map;

%% 3- Quantized Range Samples (index)
% Quantized range image (Depends on the range resolution)
RQ_idx = round(R1_map./SysSpec.dr);
RQ_idx = RQ_idx - min(min(RQ_idx))+1;
RQ_idx = bridge_mask .* RQ_idx;

%% 4- RCS for BG and CR:
% Background's maximum RCS per pixel
% RCS_BG    = 50; % e.g. 15 m2; background's RCS
RCS_BG    = .99*4*pi*((res_x*res_z)^2)/(lambda^2); %metal bridge
disp(['Background maximum RCS [each component]: ', num2str(RCS_BG)])
% CR's RCS [Trihederal]
cr_length = 0.2; % e.g. 15 cm
RCS_CR    = 12*pi*(cr_length^4)/(lambda^2); %m2
% RCS_CR    = 0;
% RCS_CR    = 500;
disp(['CR RCS [CR length=',num2str(cr_length),']: ', num2str(RCS_CR)])

%% 5- Deploy Corner reflector [in a random location]:
[SNR_noCR, slc_noCR] = SLC_scene_simulator(SysSpec, TarLoc, RCS_BG, 0, [1,1], noise_R,0);
SNR_noCR = SNR_noCR .* bridge_mask_slc;
SNR_tresh = mean(10*log10(SNR_noCR(bridge_mask==1)))+1*std(10*log10(SNR_noCR(bridge_mask==1)));
disp(['Background SNR Threshold (dB) [lower better]: ',num2str(SNR_tresh)])
RQ_idx_uniq = unique(RQ_idx);
disp(['unique ranges on target: ',num2str(length(RQ_idx_uniq)) ])
num_CR = min([15,length(RQ_idx_uniq)]);
% RQ_idx_CR = randperm(max(RQ_idx_uniq),num_CR);
RQ_idx_rand = randperm(max(RQ_idx_uniq),max(RQ_idx_uniq));
RQ_idx_CR = [];
for ii = 1:length(RQ_idx_rand)
    RQ_snr_ii  = max( 10*log10(SNR_noCR(RQ_idx == RQ_idx_rand(ii))) ); % Ba
    if RQ_snr_ii<SNR_tresh
        RQ_idx_CR = [RQ_idx_CR,RQ_idx_rand(ii)];
    end
    if length(RQ_idx_CR)==num_CR
        break
    end
end

for ii = 1:length(RQ_idx_CR)
    Loc_CR(ii,:) = [0,0];
    while Loc_CR(ii,1)<3 || Loc_CR(ii,2)<3 || Loc_CR(ii,1)>size(RQ_idx,1)-3 || Loc_CR(ii,2)>size(RQ_idx,2)-3
        % Location (Calculated based on the Quantized range array):
        % CR + targets with the same range to the sensor
        %     [row_t,col_t]  = find( RQ_idx == randi( [ min(min(RQ_idx))+1,max(max(RQ_idx))-1 ],1 ) ); % Random location!!
        [row_t,col_t]  = find( RQ_idx == RQ_idx_CR(ii) ); % Random location!!
        
        % Randomizing the locations in the array!
        rtct_dummy = [row_t,col_t];
        rtct_dummy = rtct_dummy( randperm( size( rtct_dummy,1) ),: );
        row_t      = rtct_dummy(:,1); % CR row location in quantized range (identical range cells)
        col_t      = rtct_dummy(:,2); % CR col location in quantized range (identical range cells)
        %
        Loc_CR(ii,:)     = [row_t(1), col_t(1)]; % deployed CR location
    end
end
disp(['number of deployed CRs: ',num2str(size(Loc_CR,1)) ])
%% 6- Displacement values
% Background's Displacement vec (X,Y,Z)
% epochs = 10;
% for epoch_i = 0:epochs-1
%     d_vec_BG_X(epoch_i+1,:,:) = -0*lambda/16+ ( -0.3*1e-3 ) .* ones( size(z1,1),1 )*sin( linspace(0,2*pi,size(z1,2))  );
%     d_vec_BG_Y(epoch_i+1,:,:) = 0*lambda/16 + ( 0.3*1e-3 ) .* ones( size(y1,1),1 )*sin( linspace(0,2*pi,size(y1,2))  );
%     d_vec_BG_Z(epoch_i+1,:,:) = ( -0.5*1e-3 ) .* sin( linspace(0,2*pi,size(x1,1))' )*ones( 1,size(x1,2) );
% end

epochs = 10;
for epoch_i = 0:epochs-1
    d_vec_BG_X(epoch_i+1,:,:) = -0*lambda/16+ ( -0.4*1e-3 ) .* ones( size(z1,1),1 )*sin( linspace(0,2*pi,size(z1,2))  );
    d_vec_BG_Y(epoch_i+1,:,:) = 0*lambda/16 + ( 0.4*1e-3 ) .* ones( size(y1,1),1 )*sin( linspace(0,2*pi,size(y1,2))  );
    d_vec_BG_Z(epoch_i+1,:,:) = ( -0.5*1e-3 ) .* ( linspace(0,size(x1,1),size(x1,1))'/size(x1,1) )*ones( 1,size(x1,2) );
end

% for epoch_i = 1:epochs
% d_vec_BG_X(epoch_i,:,:) = 1*lambda/8;
% d_vec_BG_Y(epoch_i,:,:) = 1*lambda/8;
% d_vec_BG_Z(epoch_i,:,:) = -1*lambda/8;
% end

% Displacement vectors [sample x epochs]
cumsumX=1e3*cumsum( squeeze(d_vec_BG_X(:,1,:)),1 )';
cumsumY=1e3*cumsum( squeeze(d_vec_BG_Y(:,1,:)),1 )';
cumsumZ=1e3*cumsum( squeeze(d_vec_BG_Z(:,:,1)),1 )';

bridge_mask = reshape(bridge_mask,1,size(bridge_mask,1),size(bridge_mask,2));
d_vec_BG_X = bridge_mask .* d_vec_BG_X;
d_vec_BG_Y = bridge_mask .* d_vec_BG_Y;
d_vec_BG_Z = bridge_mask .* d_vec_BG_Z;

d_vec_BG_X_total = 1000*squeeze(sum(d_vec_BG_X,1));d_vec_BG_X_total(bridge_mask==0)=NaN;
d_vec_BG_Y_total = 1000*squeeze(sum(d_vec_BG_Y,1));d_vec_BG_Y_total(bridge_mask==0)=NaN;
d_vec_BG_Z_total = 1000*squeeze(sum(d_vec_BG_Z,1));d_vec_BG_Z_total(bridge_mask==0)=NaN;
min_defo = min([d_vec_BG_Y_total(:)',d_vec_BG_Y_total(:)',d_vec_BG_Z_total(:)']);
max_defo = max([d_vec_BG_Y_total(:)',d_vec_BG_Y_total(:)',d_vec_BG_Z_total(:)']);

figure('Position', [0 0 1500 300]); %sgtitle('Displacement functions')
cm = colormap(jet(10)); 
for ii=1:epochs
ax(1)=subplot(1,3,1); plot( -tar_BG_X/2:res_x:tar_BG_X/2-res_x,cumsumX(:,ii),'Color',cm(ii,:) );hold on;ylim([min_defo,max_defo]);xlim([-1,1]*tar_BG_X/2);ylabel('Displacement (mm)');xlabel('X (m)')
ax(2)=subplot(1,3,2); plot( -tar_BG_X/2:res_x:tar_BG_X/2-res_x,cumsumY(:,ii),'Color',cm(ii,:) );hold on;ylim([min_defo,max_defo]);xlim([-1,1]*tar_BG_X/2);ylabel('Displacement (mm)');xlabel('X (m)')
ax(3)=subplot(1,3,3); plot( -tar_BG_Z/2:res_z:tar_BG_Z/2-res_z,cumsumZ(:,ii),'Color',cm(ii,:) );hold on;ylim([min_defo,max_defo]);xlim([-1,1]*tar_BG_Z/2);ylabel('Displacement (mm)');xlabel('Z (m)')
end
cbar = colorbar('Position',[0.93 0.168 0.022 0.7]);caxis([1,10])
cbar.Title.String = "epochs";
colormap(cbar,'jet')
if export_flag
    print(gcf,[export_directory 'Displacement functions.jpg'],'-djpeg','-r400'); 
end

figure('Position', [0 0 1500 300]);colormap( [1 1 1; jet(256)] )
subplot(1,3,1);imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1], flipud(d_vec_BG_X_total));caxis([min_defo max_defo]);axis('tight'); axis('equal');xlabel('X (m)');ylabel('Z (m)');set(gca,'YDir','normal')
subplot(1,3,2);imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1], flipud(d_vec_BG_Y_total));caxis([min_defo max_defo]);axis('tight'); axis('equal');xlabel('X (m)');ylabel('Z (m)');set(gca,'YDir','normal')
subplot(1,3,3);imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1], flipud(d_vec_BG_Z_total));caxis([min_defo max_defo]);axis('tight'); axis('equal');xlabel('X (m)');ylabel('Z (m)');set(gca,'YDir','normal')
cbar = colorbar('Position',[0.93 0.168 0.022 0.7]); 
cbar.Title.String = "mm"; 
colormap(cbar,'jet')

if export_flag
   print(gcf,[export_directory 'Reference displacement map.jpg'],'-djpeg','-r400'); 
end
%% 7- Targets' location [After Displacement]
% Background Targets
x_ts(1,:,:) = x1;
y_ts(1,:,:) = y1;
z_ts(1,:,:) = z1;
x_ts = [x_ts; x_ts+cumsum(d_vec_BG_X,1)];
y_ts = [y_ts; y_ts+cumsum(d_vec_BG_Y,1)];
z_ts = [z_ts; z_ts+cumsum(d_vec_BG_Z,1)];

% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
%% SAR & RAR SIGNALS:
%% 1- SLC
for epoch_i=1:epochs+1
    TarLoc_i =  struct;
    TarLoc_i.x = squeeze( x_ts(epoch_i,:,:) );
    TarLoc_i.y = squeeze( y_ts(epoch_i,:,:) );
    TarLoc_i.z = squeeze( z_ts(epoch_i,:,:) );
    R_map_ts(epoch_i,:,:) = sqrt(x_ts(epoch_i,:,:).^2+y_ts(epoch_i,:,:).^2+z_ts(epoch_i,:,:).^2);
    if aux_figs_flag==1 && epoch_i==1
        plot_results=1;
    else
        plot_results=0;
    end
    [SNR_ts(epoch_i,:,:),slc_ts(epoch_i,:,:)] = SLC_scene_simulator(SysSpec, TarLoc_i, RCS_BG, RCS_CR, Loc_CR, noise_R,plot_results);
end
bridge_mask_slc = reshape(bridge_mask_slc,1,size(bridge_mask_slc,1),size(bridge_mask_slc,2));
slc_ts = bridge_mask_slc .* slc_ts;
snr = 0;
slc_ts = awgn(slc_ts,snr); % add Gaussian noise with 10dB SNR


% % FROM RAW SIGNAL: VERY SLOW AND COMPUTATIONALY EXPENSIVE!!!!
% signal_SAR1 = 0;
% for ii = 1:2%size(x1,1)
%     for jj = 1:2%size(x1,2)
%         R_tar1_i = sqrt( (x1(ii,jj)-x_sensor(:)).^2 + (y1(ii,jj)-y_sensor(:)).^2 + (z1(ii,jj)-z_sensor(:)).^2 ); % Before displacement (1st epoch)
%         Amp1_i = SNR1(ii,jj);
%         signal_SAR1 = signal_SAR1 + signal_model_TS(Amp1_i, tm, fs, lambda, bw, R_tar1_i, snr);
%     end
% end

%% 2- PS Detection (Finding strong scattering points!)
slc_dB = squeeze( 10*log10(abs(slc_ts(1,:,:))) );

% [ps_r, ps_c] = find(slc_dB>0.8*max(max(slc_dB)));
% Loc_PS = [ps_r, ps_c];
% num_PS = length(ps_r);
row_PS = {};
col_PS = {};
RQ_idx_PS = 1:max(RQ_idx(:));
for ii = 1:max(RQ_idx(:))
    idx_ii = find( RQ_idx == ii );
    [row_PS{ii},col_PS{ii}] = find( RQ_idx == ii );
    [ps_r(ii,1), ps_c(ii,1)] = find( slc_dB == max(slc_dB(idx_ii)) );
    
    rcs_PS(ii,1) = slc_ts( 1,ps_r(ii), ps_c(ii) );
    rcs_clutter(ii,1) = 0;
    for jj = 1:length(row_PS{ii})
        rcs_clutter(ii,1) = rcs_clutter(ii,1)+slc_ts(1,row_PS{ii}(jj),col_PS{ii}(jj));
    end
end
Loc_PS = [ps_r, ps_c];
rcs_clutter = rcs_clutter-rcs_PS;
scr_PS = 10*log10( abs(rcs_PS)./abs(rcs_clutter) );
% if aux_figs_flag
%     figure; plot( 1:num_PS,10*log10(abs(rcs_PS)) ); hold on; plot( 1:num_PS,scr_PS ); legend('RCS','SCR')
% end

scr_thresh = 15; %dB
high_scr_idx = find( scr_PS>scr_thresh ); % Keeping PS point with SCR>5 dB
scr_PS = scr_PS(high_scr_idx);
rcs_PS = rcs_PS(high_scr_idx);
RQ_idx_PS = RQ_idx_PS(high_scr_idx);
row_PS = row_PS(high_scr_idx); 
col_PS = col_PS(high_scr_idx);
Loc_PS = Loc_PS(high_scr_idx,:);

high_snr_idx = find(10*log10(abs(rcs_PS))>SNR_tresh); % Keeping PS point with high SNR
scr_PS = scr_PS(high_snr_idx);
rcs_PS = rcs_PS(high_snr_idx);
RQ_idx_PS = RQ_idx_PS(high_snr_idx);
row_PS = row_PS(high_snr_idx); col_PS = col_PS(high_snr_idx);
Loc_PS = Loc_PS(high_snr_idx,:);
num_PS = length(Loc_PS);

disp(['Number of Detected PS points: ',num2str(num_PS)])
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
%% LOS INTERFEROMETRY:
%% 1- InSAR [ LOS Displacement (mm) ]
% dR_map_ref   = 10^3*(R2_map-R1_map); % mm
dR_map_ref   = 10^3*( diff(R_map_ts,1) ); % mm
dR_map_ref   = sum(dR_map_ref,1);
slc_ts_angle = angle(slc_ts);
dR_map_insar = -(10^3)*wrapping_operator( diff(slc_ts_angle,1) )*lambda/4/pi; % mm
dR_map_insar = bridge_mask.*sum(dR_map_insar,1);
displacement_error_los_insar_map = dR_map_ref-dR_map_insar;
disp([ "InSAR map mean displacement error (mm): ", num2str( mean(displacement_error_los_insar_map,[1,2,3]) ) ] )
if aux_figs_flag
    figure;sgtitle('LOS')
    subplot(1,2,1);imagesc(squeeze(dR_map_insar));subplot(1,2,2);imagesc(squeeze(dR_map_ref))
end
% *
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
%% 1- CR's bins in RC signals for each epoch [Range Compressed (RAR) observations of each sensor's location on SAR trajectory]
% It might take a few minutes....
x_sensor  = x_sensor(:);
z_sensor  = z_sensor(:);
for epoch_i = 1:epochs+1
    for jj = 1:num_PS
        row_t_i = row_PS{jj};
        col_t_i = col_PS{jj};
        
        RC_sig_PS_ts_i = 0;
        % Adding scattering contriutions of targets that are located in a same range bin in RAR range compressed signal:
        for ii = 1:size(row_t_i,1)
            xt_i = x_ts( epoch_i,row_t_i(ii),col_t_i(ii) );
            yt_i = y_ts( epoch_i,row_t_i(ii),col_t_i(ii) );
            zt_i = z_ts( epoch_i,row_t_i(ii),col_t_i(ii) );
            % Range between the PS and antennas' location at each step:
            R_i_ts = sqrt( (xt_i-x_sensor).^2 + (yt_i-y_sensor).^2 + (zt_i-z_sensor).^2 ); 
            % Range compressed signal:
            RC_sig_PS_ts_i = RC_sig_PS_ts_i + SNR_ts( epoch_i,row_t_i(ii),col_t_i(ii) ) .* exp( 1i*4*pi*(R_i_ts)/lambda );
        end
%         R_PS = sqrt( (x_ts(epoch_i,Loc_PS(jj,1),Loc_PS(jj,2))-x_sensor).^2 + ...
%             (y_ts(epoch_i,Loc_PS(jj,1),Loc_PS(jj,2))-y_sensor).^2 + ...
%             (z_ts(epoch_i,Loc_PS(jj,1),Loc_PS(jj,2))-z_sensor).^2 );
        
        % Complex (Real and Imaginary):
        RC_sig_PS_ts(epoch_i,jj,:) = RC_sig_PS_ts_i; % 3d Complex array: [time_samples x num_PS x radar_steps] 
    end
end

% -------------------------------------------------------------------------
%% 4- 3D Displacement [Cartesian]
SAR3dLoc(:,:,1)=SysSpec.x_n;SAR3dLoc(:,:,2)=SysSpec.y_n;SAR3dLoc(:,:,3)=SysSpec.z_n; %Radar's 3D position during data aquisition for SAR imaging. [3d array NxMx3 ]
for ps_i = 1:num_PS
    Tar3dLoc(ps_i,:) = [TarLoc_hat.x(Loc_PS(ps_i,1),Loc_PS(ps_i,2)),...
                        TarLoc_hat.y(Loc_PS(ps_i,1),Loc_PS(ps_i,2)),...
                        TarLoc_hat.z(Loc_PS(ps_i,1),Loc_PS(ps_i,2))];
end
[d_hat_3D_cartesian_dif,d_hat_3D_cartesian_ts,d_hat_3D_cartesian_total] = Displacement_vec_Cartesian(Tar3dLoc,SAR3dLoc,RC_sig_PS_ts,lambda,0);
% OUTPUTS:
% d_hat_3D_cartesian_dif: amount of targets' differential displacemet vector at each time (epoch). 3d array: [epochs x targets x xyz]
% d_hat_3D_cartesian_ts: cumulative displacemet vector at each time.                               3d array: [epochs x targets x xyz]
% d_hat_3D_cartesian_total: total displacemet vector after time-series.                            2d array: [targets x xyz]

% REFERENCE:
for epoch_i = 1:epochs
    for ps_i = 1:num_PS
        d_vec_PS_dif(epoch_i,ps_i,:) = [d_vec_BG_X(epoch_i,Loc_PS(ps_i,1),Loc_PS(ps_i,2)); ...
                                    d_vec_BG_Y(epoch_i,Loc_PS(ps_i,1),Loc_PS(ps_i,2)); ...
                                    d_vec_BG_Z(epoch_i,Loc_PS(ps_i,1),Loc_PS(ps_i,2))];
    end
end
d_vec_PS_ts    = cumsum(d_vec_PS_dif,1);
d_vec_PS_total = squeeze( d_vec_PS_ts(end,:,:) ); % targets x xyz
if num_PS==1
    d_vec_PS_total = d_vec_PS_total';
end

% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------------------------


%% RESULTS:

%%  METRICS
error_d_hat_3D_cartesian_total = d_hat_3D_cartesian_total-d_vec_PS_total; % targets x xyz
norm_d_vec_total = sum(d_vec_PS_total.^2,2).^0.5; % targets x 1
norm_estimated_d_vec_total = sum(d_hat_3D_cartesian_total.^2,2).^0.5; % targets x 1
error_norm_d_vec_total = norm_d_vec_total-norm_estimated_d_vec_total; % targets x 1
rmse_abs_defo_error_total = sqrt(mean((error_norm_d_vec_total).^2)); %rmse abs (m) 1x1
rmse_vec_total = sqrt( mean(error_d_hat_3D_cartesian_total.^2,1) )'; %rmse x,y,z (m) xyz x 1
disp([" 3D displacement's RMSE in mm (dx,dy,dz): "; num2str(1000*rmse_vec_total)])
error_angle_d_vec_total = acosd( sum(d_vec_PS_total.*d_hat_3D_cartesian_total,2)./(norm_estimated_d_vec_total.*norm_d_vec_total) ); % targets x 1
% error_angle_d_vec_total = acosd( dot(d_vec_PS_total'./norm_d_vec_total', d_hat_3D_cartesian_total'./norm_estimated_d_vec_total') )';
rmse_angle_d_vec_total = sqrt(mean(error_angle_d_vec_total.^2)); % 1x1


error_d_hat_3D_cartesian_dif = d_hat_3D_cartesian_dif-d_vec_PS_dif; % epochs x targets x xyz
norm_d_vec_dif = sum(d_vec_PS_dif.^2,3).^0.5; % epochs x targets
norm_estimated_d_vec_dif = sum(d_hat_3D_cartesian_dif.^2,3).^0.5; % epochs x targets
error_norm_d_vec_dif = norm_d_vec_dif-norm_estimated_d_vec_dif; % epochs x targets
rmse_abs_defo_error_dif = sqrt(mean((error_norm_d_vec_dif).^2,2)); %rmse abs (m) epochs x 1
rmse_vec_dif = squeeze( sqrt(mean(error_d_hat_3D_cartesian_dif.^2,2)) ); %rmse x,y,z (m) xyz x epochs
if epochs>1
    rmse_vec_dif = rmse_vec_dif'; %rmse x,y,z (m) xyz x epochs
end
disp([" 3D displacement's RMSE in mm (dx,dy,dz): "; num2str(1000*rmse_vec_dif)])
for epoch_i = 1:epochs
    d_vec_PS_dif_i = squeeze(d_vec_PS_dif(epoch_i,:,:))';
    d_hat_3D_cartesian_dif_i = squeeze(d_hat_3D_cartesian_dif(epoch_i,:,:))';
    error_angle_d_vec_dif(epoch_i,:) = acosd( dot(d_vec_PS_dif_i./norm_d_vec_dif(epoch_i,:), d_hat_3D_cartesian_dif_i./norm_estimated_d_vec_dif(epoch_i,:)) )'; % epochs x targets
end
rmse_angle_d_vec_dif = sqrt(mean(error_angle_d_vec_dif.^2,2)); % epochs x 1

%%
%% = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
%% PLOTS
    slc1 = squeeze(slc_ts(1,:,:));
    slc2 = squeeze(slc_ts(2,:,:));
    cbar_min = 0;
    cbar_max = 10*log10(max(max(abs(slc1))));
% % % if aux_figs_flag
% % %     %% 0- General overview of scene simulations
% % %     figure;set(gca,'YDir','normal')
% % %     
% % %     ax1 = subplot(4,8,[1:4]);   imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud(10*log10(abs(slc1))), [ cbar_min, cbar_max ] ); 
% % %     ax = gca;ax.YDir= 'normal'; axis('equal'); title('slc-abs (dB)'); xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
% % %     ax2 = subplot(4,8,[9:12]);  imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud(angle(slc1)));axis('equal'); title('slc-angle');
% % %     ax = gca;ax.YDir= 'normal'; xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
% % %     ax3 = subplot(4,8,[17:20]); imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud((10^3)*wrapping_operator(angle(slc1)-angle(slc2))*lambda/4/pi) );
% % %     ax = gca;ax.YDir= 'normal';axis('equal'); title('LOS displacement map (mm)'); xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
% % %     
% % %     linkaxes([ax1, ax2, ax3]);
% % %     subplot(4,8,[5,6]);   hist( angle(slc1(:)),5000 );         title('Phase distribution (Rad)')
% % %     subplot(4,8,[7,8]);   hist( 10*log10(abs(slc1(:))),5000 ); title('Amplitude distribution (dB)')
% % %     subplot(4,8,[13,14]); hist( real(slc1(:)),50000 );         title('Real distribution')
% % %     subplot(4,8,[15,16]); hist( imag(slc1(:)),50000 );         title('Imaginary distribution')
% % %     subplot(4,8,[25:28]); imagesc([-tar_BG_X/2 tar_BG_X/2],[0 tar_BG_Z],flipud(RQ_idx));
% % %     ax = gca;ax.YDir= 'normal';axis('equal'); colormap('jet');title('iso-range contours');xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
% % %     subplot(4,8,[23,24,31,32]);
% % %     scatter3( x1(:),y1(:),z1(:),'filled','k')
% % %     hold on
% % %     scatter3( x_sensor(:),y_sensor*x_sensor(:),z_sensor(:), 30,'filled')
% % %     lgd=legend('Scene area', 'Sensor steps'); %axis equal
% % %     lgd.Location = 'southoutside';
% % %     title('Scene geometry'); xlabel('Cross-range (x)'); ylabel('Range (y)');zlabel('Elevation (z)');
% % %     
% % %     subplot(4,8,[21,22,29,30]); %plot( [0,real(tt)], [0,imag(tt)] );%axis('equal')
% % %     hold on
% % %     
% % %     for ii = 1:size(row_t,1)
% % %         if ii==1
% % %             plot( [0, real( slc1(row_t(ii),col_t(ii)) )], [0,imag( slc1(row_t(ii),col_t(ii)) )] )
% % %             base_real=0; base_imag=0;
% % %         else
% % %             base_real=base_real+real( slc1(row_t(ii-1),col_t(ii-1)) );
% % %             base_imag=base_imag+imag( slc1(row_t(ii-1),col_t(ii-1)) );
% % %             plot( [base_real, base_real+real( slc1(row_t(ii),col_t(ii)) )], [base_imag,base_imag+imag( slc1(row_t(ii),col_t(ii)) )] )
% % %         end
% % %     end
% % %     title("Targets' interference for CR's pixel")
% % %     xlabel('Real'); ylabel('Imaginary')
% % % end

%% Figure 1- Simulated slc
xfig        = 0;   % Screen position
yfig        = 0;   % Screen position
fontsizefig = 14; fontname    = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',12);

figure('Position', [xfig yfig 1700 1000]);set(gca,'YDir','normal')
sgtitle('Simulations')
ax1 = subplot(2,2,1); imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud(10*log10(abs(slc1))), [ cbar_min, cbar_max ] );
ax = gca;ax.YDir= 'normal'; axis('tight'); axis('equal'); title('slc-abs (dB)'); xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
ax2 = subplot(2,2,2); imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud(angle(slc1)));
axis('tight');axis('equal'); title('slc-angle'); ax = gca;ax.YDir= 'normal'; xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
ax3 = subplot(2,2,3); imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud((10^3)*wrapToPi(angle(slc1)-angle(slc2))*lambda/4/pi) );
ax = gca;ax.YDir= 'normal'; axis('tight');axis('equal'); title('LOS displacement map (mm)'); xlabel('Azimuth (x)'); ylabel('Elevation (z)'); colorbar()
% subplot(2,2,4); imagesc(tar_BG_X/2*[-1 1],tar_BG_Z/2*[-1 1],flipud(RQ_idx));
subplot(2,2,4); RQ_idx(RQ_idx==0)=NaN;contourf(flipud(RQ_idx),max(max(RQ_idx)));
colormap('jet'); title('iso-range contours');xlabel('Azimuth (x)'); ylabel('Elevation (z)'); ax = gca;ax.YDir= 'normal'; % axis('auto');  %colorbar()
if export_flag
   print(gcf,[export_directory 'SLC simulations.jpg'],'-djpeg','-r600'); 
end
   print(gcf,[export_directory 'Scene simulation process.jpg'],'-djpeg','-r600'); 

%% Figure 2- Histogram
% % % Displacement vectors
figure('Position', [0 0 1200 1200]);
sgtitle('Histograms of displacement errors')

xlim_hist_max = []; xlim_hist_min = [];
ylim_hist_max = []; ylim_hist_min = [];
ylim_hist_angle_max = []; ylim_hist_angle_min = [];
for epoch_i = 1:epochs
    subplot( epochs,5,(epoch_i-1)*5+1 )
    hist_dif(epoch_i,1,:) = histogram( 1000*error_d_hat_3D_cartesian_dif(epoch_i,:,1), round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );
    
    subplot( epochs,5,(epoch_i-1)*5+2 )
    hist_dif(epoch_i,2,:) = histogram( 1000*error_d_hat_3D_cartesian_dif(epoch_i,:,2), round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );
    subplot( epochs,5,(epoch_i-1)*5+3 )
    hist_dif(epoch_i,3,:) = histogram( 1000*error_d_hat_3D_cartesian_dif(epoch_i,:,3), round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );
    subplot( epochs,5,(epoch_i-1)*5+4 )
    hist_dif(epoch_i,4,:) = histogram( 1000*(error_norm_d_vec_dif(epoch_i,:)),         round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );
    subplot( epochs,5,(epoch_i-1)*5+5 )
    hist_dif(epoch_i,5,:) = histogram( error_angle_d_vec_dif(epoch_i,:),               round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );
    ylim_hist_max = [ylim_hist_max,max(hist_dif(epoch_i,1,:).Values)];
    ylim_hist_max = [ylim_hist_max,max(hist_dif(epoch_i,2,:).Values)];
    ylim_hist_max = [ylim_hist_max,max(hist_dif(epoch_i,3,:).Values)];
    ylim_hist_max = [ylim_hist_max,max(hist_dif(epoch_i,4,:).Values)];
    ylim_hist_angle_max = [ylim_hist_angle_max,max(hist_dif(epoch_i,5,:).Values)];
end
xlim_defo0 = 1100*max([max(abs(error_d_hat_3D_cartesian_dif(:))),abs(error_norm_d_vec_dif(:))',0*1e-4]); xlim_defo = xlim_defo0*[ -1 , 1 ];
xlim_defo_angle0 = 1.1*max([0*5,abs(error_angle_d_vec_dif(:))']); xlim_defo_angle = xlim_defo_angle0*[ -1 , 1 ];
ylim_hist = [ 0 , 1.1*max(ylim_hist_max) ];
ylim_hist_angle = [ 0 , 1.1*max(ylim_hist_angle_max) ];
for epoch_i = 1:epochs
    subplot( epochs,5,(epoch_i-1)*5+1 )
    xlim( xlim_defo );ylim( ylim_hist );ylabel('Frequency');set(gca,'YGrid','on')
    subplot( epochs,5,(epoch_i-1)*5+2 )
    xlim( xlim_defo );ylim( ylim_hist );set(gca,'YGrid','on')
    subplot( epochs,5,(epoch_i-1)*5+3 )
    xlim( xlim_defo );ylim( ylim_hist );set(gca,'YGrid','on')
    subplot( epochs,5,(epoch_i-1)*5+4 )
    xlim( xlim_defo );ylim( ylim_hist );set(gca,'YGrid','on')
    subplot( epochs,5,(epoch_i-1)*5+5 )
    xlim( xlim_defo_angle );ylim( ylim_hist_angle );set(gca,'YGrid','on')
end
subplot( epochs,5,1 );title('X error distribution (mm)')
subplot( epochs,5,2 );title('Y error distribution (mm)')
subplot( epochs,5,3 );title('Z error distribution (mm)')
subplot( epochs,5,4 );title('Absolute error distribution (mm)')
subplot( epochs,5,5 );title('angle error distribution (deg)')

subplot( epochs,5,(epochs-1)*5+1 );xlabel('Error distribution (mm)')
subplot( epochs,5,(epochs-1)*5+2 );xlabel('Error distribution (mm)')
subplot( epochs,5,(epochs-1)*5+3 );xlabel('Error distribution (mm)')
subplot( epochs,5,(epochs-1)*5+4 );xlabel('Error distribution (mm)')
subplot( epochs,5,(epochs-1)*5+5 );xlabel('Error distribution (deg)')
if export_flag
   print(gcf,[export_directory 'Histograms of displacement errors for each TS.jpg'],'-djpeg','-r600'); 
end
%% Figure 3 
figure('Position', [0 0 1200 1200]);
sgtitle('displacement errors in TS')
xax = linspace(1,epochs,epochs);
subplot(3,2,1); plot(1000*rmse_vec_dif(1,:),'ko-','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',5);xticks(xax)
title('X displacement error in time-series');xlabel('Time samples');ylabel('mm')
subplot(3,2,3); plot(1000*rmse_vec_dif(2,:),'ko-','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',5);xticks(xax)
title('Y displacement error in time-series');xlabel('Time samples');ylabel('mm')
subplot(3,2,5); plot(1000*rmse_vec_dif(3,:),'ko-','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',5);xticks(xax)
title('Z displacement error in time-series');xlabel('Time samples');ylabel('mm')
subplot(3,2,2); plot(1000*rmse_abs_defo_error_dif,'ko-','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',5);xticks(xax)
title('Absolute displacement error in time-series');xlabel('Time samples');ylabel('mm')
subplot(3,2,4); plot(rmse_angle_d_vec_dif,'ko-','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2,'MarkerSize',5);xticks(xax)
title('Angular displacement error in time-series');xlabel('Time samples');ylabel('deg')
if export_flag
   print(gcf,[export_directory 'plots of displacement errors.jpg'],'-djpeg','-r600'); 
end

%% Figure 4
figure('Position', [0 0 1200 1200]);
sgtitle('Histograms of displacement errors (all samples together)')
subplot(3,2,1)
hist1 = histogram( 1000*error_d_hat_3D_cartesian_dif(:,:,1), round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );title('X error distribution (mm)');
subplot(3,2,3)
hist2 = histogram( 1000*error_d_hat_3D_cartesian_dif(:,:,2), round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );title('Y error distribution (mm)')
subplot(3,2,5)
hist3 = histogram( 1000*error_d_hat_3D_cartesian_dif(:,:,3), round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );title('Z error distribution (mm)')
subplot(3,2,2)
hist4 = histogram( 1000*(error_norm_d_vec_dif),              round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );title('abs error distribution (mm)')
subplot(3,2,4)
hist5 = histogram( error_angle_d_vec_dif,                    round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',2 );title('angle error distribution (deg)')

xlim_defo0 = 1100*max([max(abs(error_d_hat_3D_cartesian_dif(:))),abs(error_norm_d_vec_dif(:))']+0*.1e-4); xlim_defo = xlim_defo0*[ -1 , 1 ];
xlim_defo_angle0 = 1.1*max(0+abs(error_angle_d_vec_dif(:))); xlim_defo_angle = xlim_defo_angle0*[ -1 , 1 ];
ylim_hist = [ 0 , 1.1*max([hist1.Values,hist2.Values,hist3.Values,hist4.Values,hist5.Values]) ];
ylim_hist_angle = [ 0 , 1.1*max([hist1.Values,hist2.Values,hist3.Values,hist4.Values,hist5.Values]) ];

subplot(3,2,1);xlim( xlim_defo );ylim( ylim_hist );ylabel('Frequency');xlabel('Error distribution (mm)');set(gca,'YGrid','on')
subplot(3,2,2);xlim( xlim_defo );ylim( ylim_hist );ylabel('Frequency');xlabel('Error distribution (mm)');set(gca,'YGrid','on')
subplot(3,2,3);xlim( xlim_defo );ylim( ylim_hist );ylabel('Frequency');xlabel('Error distribution (mm)');set(gca,'YGrid','on')
subplot(3,2,5);xlim( xlim_defo );ylim( ylim_hist );ylabel('Frequency');xlabel('Error distribution (mm)');set(gca,'YGrid','on')
subplot(3,2,4);xlim( xlim_defo_angle );ylim( ylim_hist_angle );ylabel('Frequency');xlabel('Error distribution (deg)');set(gca,'YGrid','on')
if export_flag
   print(gcf,[export_directory 'Histograms of displacement errors for each TS all samples together.jpg'],'-djpeg','-r600'); 
end

%% Figure 5
figure('Position', [0 0 1200 1200]);
sgtitle('Histograms of cumulative displacement errors')
subplot(3,2,1)
hist1 = histogram( 1000*error_d_hat_3D_cartesian_total(:,1),    round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',1 );title('X error distribution (mm)')
subplot(3,2,3)
hist2 = histogram( 1000*error_d_hat_3D_cartesian_total(:,2),    round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',1 );title('Y error distribution (mm)')
subplot(3,2,5)
hist3 = histogram( 1000*error_d_hat_3D_cartesian_total(:,3),    round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',1 );title('Z error distribution (mm)')
subplot(3,2,2)
hist4 = histogram( 1000*(error_norm_d_vec_total),               round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',1 );title('abs error distribution (mm)')
% [f,xi] = ksdensity(1000*(norm_d_vec-norm_estimated_d_vec));
% hold on;plot(xi,f)
subplot(3,2,4)
hist5 = histogram( error_angle_d_vec_total,                     round(num_PS/.1),'EdgeColor', [0.53,0.81,1], 'FaceColor', [0.53,0.81,1], 'FaceAlpha', 1, 'LineWidth',1 );title('angle error distribution (deg)')


xlim_defo0 = 1100*max([max(abs(error_d_hat_3D_cartesian_total(:))),abs(error_norm_d_vec_total(:))',0*1e-4]); xlim_defo = xlim_defo0*[ -1 , 1 ];
xlim_defo_angle0 = 1.1*max([abs(error_angle_d_vec_total(:))']); xlim_defo_angle = xlim_defo_angle0*[ -1 , 1 ];
ylim_hist = [ 0 , 1.1*max([hist1.Values,hist2.Values,hist3.Values,hist4.Values,hist5.Values]) ];
ylim_hist_angle = [ 0 , 1.1*max([hist1.Values,hist2.Values,hist3.Values,hist4.Values,hist5.Values]) ];

subplot(3,2,1);xlim( xlim_defo );ylim( ylim_hist )
subplot(3,2,2);xlim( xlim_defo );ylim( ylim_hist )
subplot(3,2,3);xlim( xlim_defo );ylim( ylim_hist )
subplot(3,2,5);xlim( xlim_defo );ylim( ylim_hist )
subplot(3,2,4);xlim( xlim_defo_angle );ylim( ylim_hist_angle )
if export_flag
   print(gcf,[export_directory 'Histograms of cumulative displacement errors.jpg'],'-djpeg','-r600'); 
end

%% Figure 6-8 - Time series displacements
% figure('Position', [0 0 1200 1200]);
% sgtitle('Time series displacements')
% ax_ts = 1:1:epochs;
% ax_ts2 = [ax_ts,fliplr(ax_ts)];
% for PS_i = 1:num_PS
%     for xyz = 1:3  %xyz directions
%         subplot(num_PS,3,(PS_i-1)*3+xyz); curve_dev = squeeze(d_hat_3D_cartesian_ts(:,PS_i,xyz))+std(squeeze(d_hat_3D_cartesian_ts(:,PS_i,xyz)))*[1 -1];
%         inBetween = [curve_dev(:,1)', fliplr(curve_dev(:,2)')];
%         fill(ax_ts2, inBetween,[.7,.7,.7],'FaceAlpha',.5,'EdgeColor', 'w'); hold on
%         scatter(ax_ts,d_vec_PS_ts(:,PS_i,xyz),10,'ko','Filled');
%         plot( d_hat_3D_cartesian_ts(:,PS_i,xyz),'k' )
%     end
% end

XYZ = ['X','Y','Z'];
for xyz = 1:3
    figure('Position', [0 0 1200 1300]);
    sgtitle(['Time series displacements in ', XYZ(xyz)])
    ax_ts = 1:1:epochs;
    ax_ts2 = [ax_ts,fliplr(ax_ts)];
    subp_col = floor(sqrt(num_PS));
    subp_row = ceil(num_PS/subp_col)+1;
    for PS_i = 1:num_PS
        subplot(subp_row,subp_col,PS_i);
        curve_dev = 1000*squeeze(d_hat_3D_cartesian_ts(:,PS_i,xyz))+1000*std(squeeze(d_hat_3D_cartesian_ts(:,PS_i,xyz)))*[1 -1];
        inBetween = [curve_dev(:,1)', fliplr(curve_dev(:,2)')];
        fill(ax_ts2, inBetween,[.7,.7,.7],'FaceAlpha',.5,'EdgeColor', 'w'); hold on
        plot( 1000*d_hat_3D_cartesian_ts(:,PS_i,xyz),'k' )
        scatter(ax_ts,1000*d_vec_PS_ts(:,PS_i,xyz),'kx');
    end
    
    % add legend
    Lgnd = legend({'Measured displacement buffer','Measured displacement','Actual displacement'},'FontSize',14);
    Lgnd.Position(1) = 0.65;
    Lgnd.Position(2) = 0.1;
    if export_flag
        print(gcf,[export_directory, 'Time series results for ',XYZ(xyz),'.jpg'],'-djpeg','-r600');
    end

end

%% Figure 9- Abstract plot (Only for the method 2:Cartesian)
figure('Position', [0 0 1700 1000]);
sgtitle('Cumulative TS results')
lim = 1.2*max(max(1000*abs(error_d_hat_3D_cartesian_total)));
subplot(2,2,1);bar(categorical({'x','y','z','abs'}),1000*[rmse_vec_total; rmse_abs_defo_error_total],0.4,'k'); ylim([-lim lim]);
title(['Absolute RMSE: Overall_{(mm)}=',num2str(1000*rmse_abs_defo_error_total)]);xlabel('direction');ylabel('displacement (mm)');set(gca,'Fontsize',fontsizefig)
% %
subplot(2,2,2);
fimplicit( @(x,y) x.^2 + y.^2 -1, 'k', 'LineWidth',4 ); axis equal; axis('off');xlim([-1.3,1.3]);ylim([-1.3,1.3])
hold on; plot( [0,0],[0,1], '--k', 'LineWidth',2 ); hold on; plot( [0,0],[0,-1], '--k', 'LineWidth',.8 ); hold on; plot( [-1,1], [0,0],'--k', 'LineWidth',.8 )
% title(['Displacement error angle(deg): ', num2str(defo_error_angle)]);
hold on;
text( 0, 1.1, '0', 'FontSize', 12, 'Color', 'k' );
text( 1.1, 0, '90', 'FontSize', 12, 'Color', 'k' );
text( -0.1, -1.1, '180', 'FontSize', 12, 'Color', 'k' );
text( -1.3, 0, '270', 'FontSize', 12, 'Color', 'k' );
for ii = 1:num_PS
    plot([0,sind(error_angle_d_vec_total(ii))] , [0,cosd(error_angle_d_vec_total(ii))],'--r', 'LineWidth',2 )
    hold on; %text( (1.1)*sind(defo_error_angle), (1.1)*cosd(defo_error_angle), num2str(defo_error_angle), 'FontSize', 12, 'Color', 'r' );
end
title( ['Cumulative angle error | RMSE (deg): ', num2str(rmse_angle_d_vec_total)] )
% %
subplot(2,2,3)
defo_dsp=quiver3(zeros(size(d_vec_PS_total,1),1),zeros(size(d_vec_PS_total,1),1),zeros(size(d_vec_PS_total,1),1), ...
    1000*d_vec_PS_total(:,1),1000*d_vec_PS_total(:,2),1000*d_vec_PS_total(:,3), 'b','LineWidth',1, 'ShowArrowHead', 'on');
defo_dsp.MaxHeadSize=0.5;
hold on
defo_dsp_hat2=quiver3(zeros(size(d_hat_3D_cartesian_total,1),1),zeros(size(d_hat_3D_cartesian_total,1),1),zeros(size(d_hat_3D_cartesian_total,1),1), ...
    1000*d_hat_3D_cartesian_total(:,1),1000*d_hat_3D_cartesian_total(:,2),1000*d_hat_3D_cartesian_total(:,3), 'g','LineWidth',1, 'ShowArrowHead', 'on');
defo_dsp_hat2.MaxHeadSize=0.5;
axis equal
x_lim=1.1*max(1000*abs([d_vec_PS_total(:);d_hat_3D_cartesian_total(:)]));
y_lim=1.1*max(1000*abs([d_vec_PS_total(:);d_hat_3D_cartesian_total(:)]));
z_lim=1.1*max(1000*abs([d_vec_PS_total(:);d_hat_3D_cartesian_total(:)]));
xlim( [-x_lim x_lim] );ylim( [-y_lim y_lim] );zlim( [-z_lim z_lim] )
legend('real cumulative displacement', 'estimated cumulative displacemenet');title('Cumulative Displacement Vectors')
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)')
set(gca,'Fontsize',fontsizefig)
% %
subplot(2,2,4);
error_dsp1=quiver3(zeros(size(d_vec_PS_total,1),1),zeros(size(d_vec_PS_total,1),1),zeros(size(d_vec_PS_total,1),1), ...
    error_d_hat_3D_cartesian_total(:,1),error_d_hat_3D_cartesian_total(:,2),error_d_hat_3D_cartesian_total(:,3), 'r','LineWidth',1, 'ShowArrowHead', 'on');
error_dsp1.MaxHeadSize=0.5; %text(error_d_hat_3D_cartesian(1)+0.001,error_d_hat_3D_cartesian(2)+0.001,error_d_hat_3D_cartesian(3)+0.001, 'error');
axis equal
legend('error vector');title('Cumulative Displacement error vector')
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
lim=1.1*max(max(abs(error_d_hat_3D_cartesian_total)));
xlim( [-lim lim] );ylim( [-lim lim] );zlim( [-lim lim] )

if export_flag
   print(gcf,[export_directory 'Abstract results.jpg'],'-djpeg','-r600'); 
end

%% Figure 10- Error relations (SCR vs errors)
figure('Position', [0 0 1700 1000]);
sgtitle('Relations between SCR and estimated cumulated displacement errors')

for ii = 1:3 %xyz
    subplot(3,2,2*ii-1);
    reg_model = fitlm( 10*log10(abs(scr_PS)),1000*(abs(error_d_hat_3D_cartesian_total(:,ii))) );
    reg_param = reg_model.Coefficients.Estimate;
    plot(10*log10(abs(scr_PS)),1000*(abs(error_d_hat_3D_cartesian_total(:,ii))),'xb'); hold on
    plot(10*log10(abs(scr_PS)),reg_param(1)+10*log10(abs(scr_PS))*reg_param(2),'r');
    xlabel('SCR (dB)'); ylabel('Absolute error (mm)','FontSize',12);legend({'PS points','linear fit'},'FontSize',12)
end
subplot(3,2,1);title('SCR vs error in X axis');
subplot(3,2,3);title('SCR vs error in Y axis');
subplot(3,2,5);title('SCR vs error in Z axis');
subplot(3,2,4);
reg_model = fitlm( 10*log10(abs(scr_PS)),1000*(abs(error_norm_d_vec_total)) );
reg_param = reg_model.Coefficients.Estimate;
plot(10*log10(abs(scr_PS)),1000*(abs(error_norm_d_vec_total)),'xb'); hold on
plot(10*log10(abs(scr_PS)),reg_param(1)+10*log10(abs(scr_PS))*reg_param(2),'r');legend('PS points','linear fit')
xlabel('SCR (dB)'); ylabel('Absolute error (mm)'); title('SCR vs error in absolute displacement')
subplot(3,2,6);
reg_model = fitlm( 10*log10(abs(scr_PS)),abs(error_angle_d_vec_total) );
reg_param = reg_model.Coefficients.Estimate;
plot(10*log10(abs(scr_PS)),error_angle_d_vec_total,'xb'); hold on
plot(10*log10(abs(scr_PS)),reg_param(1)+10*log10(abs(scr_PS))*reg_param(2),'r');legend('PS points','linear fit')
xlabel('SCR (dB)'); ylabel('Absolute angle error (deg)'); title('SCR vs Angular error')

if export_flag
   print(gcf,[export_directory ' Error relations (SCR vs errors).jpg'],'-djpeg','-r600'); 
end

%% Figure 11- Error relations (RCS vs errors)
figure('Position', [0 0 1700 1000]);
sgtitle('Relations between RCS and estimated cumulated displacement errors')

for ii = 1:3 %xyz
    subplot(3,2,2*ii-1);
    reg_model = fitlm( 10*log10(abs(rcs_PS)),1000*(abs(error_d_hat_3D_cartesian_total(:,ii))) );
    reg_param = reg_model.Coefficients.Estimate;
    plot(10*log10(abs(rcs_PS)),1000*(abs(error_d_hat_3D_cartesian_total(:,ii))),'xb'); hold on
    plot(10*log10(abs(rcs_PS)),reg_param(1)+10*log10(abs(rcs_PS))*reg_param(2),'r');
    xlabel('RCS (dB)'); ylabel('Absolute error (mm)','FontSize',12);legend({'PS points','linear fit'},'FontSize',12)
end
subplot(3,2,1);title('RCS vs error in X axis');
subplot(3,2,3);title('RCS vs error in Y axis');
subplot(3,2,5);title('RCS vs error in Z axis');
subplot(3,2,4);
reg_model = fitlm( 10*log10(abs(rcs_PS)),1000*(abs(error_norm_d_vec_total)) );
reg_param = reg_model.Coefficients.Estimate;
plot(10*log10(abs(rcs_PS)),1000*(abs(error_norm_d_vec_total)),'xb'); hold on
plot(10*log10(abs(rcs_PS)),reg_param(1)+10*log10(abs(rcs_PS))*reg_param(2),'r');legend('PS points','linear fit')
xlabel('RCS (dB)'); ylabel('Absolute error (mm)'); title('RCS vs error in absolute displacement')
subplot(3,2,6);
reg_model = fitlm( 10*log10(abs(rcs_PS)),abs(error_angle_d_vec_total) );
reg_param = reg_model.Coefficients.Estimate;
plot(10*log10(abs(rcs_PS)),error_angle_d_vec_total,'xb'); hold on
plot(10*log10(abs(rcs_PS)),reg_param(1)+10*log10(abs(rcs_PS))*reg_param(2),'r');legend('PS points','linear fit')
xlabel('RCS (dB)'); ylabel('Absolute angle error (deg)'); title('RCS vs Angular error')

if export_flag
   print(gcf,[export_directory ' Error relations (RCS vs errors).jpg'],'-djpeg','-r600'); 
end


%%-------------------------------------------------------------------------
%%
%%
%% FUNCTIONS
%%
function slc_resampled = resample_slc(slc, sampling, windowing)

if windowing==0
    slc_resampled = ifft2(fft2(slc),sampling(1),sampling(2));
else
    slc_resampled = fft2(slc);
    han1 = hann(size(slc_resampled,1)); han1 = repmat(han1,1,size(slc_resampled,2));
    slc_resampled = slc_resampled.*han1;
    han2 = hann(size(slc_resampled,2))'; han2 = repmat(han2,size(slc_resampled,1),1);
    slc_resampled = slc_resampled.*han2;
    slc_resampled = ifft2(slc_resampled,sampling(1),sampling(2));
end

%% NO FFT-SHIFTS!!
% resam_slc = fftshift(fft2(slc));
% resam_slc = fftshift(ifft2( resam_slc, sampling(1),sampling(2) ));
end
