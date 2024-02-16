% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

% * This code simulates a 4D array of TS of raw cube data [r x az x el x time]
% * For further processing, including SAR and InSAR:
% * run "threeD_dsplcmnt_Part2_processing.m" after running this code!



clear;
clc;
close all
addpath('./src')
c = physconst('LightSpeed');
%% Parameters:
fc = 77e9; lambda = c/fc;
bw = 1e9;
rr = c/2/bw;
T = 60e-6; cr = bw/T;
snr = 20;
Nr = 256;          % Number of range samples
Na = 101;          % Number of azimuth samples
Ne = 11;        % Number of elevation samples
Nts= 10;           % Number of TS data
d_az = lambda/2;
d_el = 2*lambda/2;

%% Radar position:
X_rad   = d_az*[-floor(Na/2):floor(Na/2)]; %d_az*(0:Na-1)-d_az*Na/2;
Y_rad0  = zeros(Na, Ne);
Z_rad0  = d_el*[-floor(Ne/2):floor(Ne/2)]; %d_el*(0:Ne-1)-d_el*Ne/2;

% Define rotation around X axis (nadir or incidence angle):
nadir   = -0;
rot_mat = [1 0 0; 0 cosd(nadir) -sind(nadir); 0 sind(nadir) cosd(nadir)];
z_shift = 0;

% Sensor's moving surface:
Y_rad = Y_rad0*cosd(nadir)-Z_rad0*sind(nadir);         % shape: Na*Ne
Z_rad = z_shift+Y_rad0*sind(nadir)+Z_rad0*cosd(nadir); % shape: Na*Ne
X_rad = X_rad(:)*ones(1,Ne);                           % shape: Na*Ne

figure('Position', [500 500 800 500]); 
% scatter(X_rad(:),Z_rad(:)); axis equal; title("Sensor's trajectory"); xlabel('Azimuth (m)'); ylabel('Elevation (m)')
scatter3(Y_rad(:),X_rad(:),Z_rad(:)); axis auto; title("Sensor's trajectory"); xlabel('Range (m)'); ylabel('Azimuth (m)'); zlabel('Elevatio (m)')

% Reshape it to 3D cube:
X_rad= repmat( reshape(X_rad, 1,Na,Ne), Nr,1,1); % shape: Nr*Na*Ne
Y_rad= repmat( reshape(Y_rad, 1,Na,Ne), Nr,1,1); % shape: Nr*Na*Ne
Z_rad= repmat( reshape(Z_rad, 1,Na,Ne), Nr,1,1); % shape: Nr*Na*Ne

% Add stepper fluctuation function:
xrad_err=0;2*1e-5;
yrad_err=0;2*1e-5;
zrad_err=0;2*1e-5;

%% Target position:
% Polar
% R_tar=[2]; theta_tar=[20.2]; ph_tar=[-20.5];
% X_tar=R_tar.*sind(theta_tar).*cosd(ph_tar);
% Y_tar=R_tar.*cosd(theta_tar).*cosd(ph_tar);
% Z_tar=R_tar.*sind(ph_tar);

% % Cartesian:
X_tar=[ 5,  7];
Y_tar=[25, 24];
Z_tar=[ 5,  3];

% X_tar=[ 5];
% Y_tar=[25];
% Z_tar=[ 5];

%% Displacement per TS:
dX_tar = [0.5e-3,  0];
dY_tar = [-0.5e-3, 0];
dZ_tar = [0.5e-3,  0];
% dX_tar = [0.5e-3];
% dY_tar = [-0.5e-3];
% dZ_tar = [0.5e-3];

num_tar=length(X_tar);

t = linspace(0,T,Nr)'; % Pulse (sweep) time axis
t = repmat(t, 1, Na, Ne); % shape: Nr*Na*Ne

%% Raw cube SIGNAL:
for ts_ii = 1:Nts
    % Add stepper fluctuation:
    X_rad2 = X_rad + xrad_err*(rand(1, Na, Ne)-0.5);
    Y_rad2 = Y_rad + yrad_err*(rand(1, Na, Ne)-0.5);
    Z_rad2 = Z_rad + zrad_err*(rand(1, Na, Ne)-0.5);

    %Rail's starting point error (spe) in X axis (unsynchronization):
    spe_shift=1; % shift in the steps (e.g., 1 step)
    spe_X(:,ts_ii) = 0*(-1.^randi(2,Ne,1).*randi(spe_shift,Ne,1).*lambda/8); 
    X_rad2 = X_rad2+reshape(spe_X(:,ts_ii),1,1,Ne); % + Starting point error (mismatch)
    
    cube_3d_ii=0;
    for tar_ii = 1:num_tar
        X_tarii = X_tar(tar_ii)+(ts_ii-1)*dX_tar(tar_ii);
        Y_tarii = Y_tar(tar_ii)+(ts_ii-1)*dY_tar(tar_ii);
        Z_tarii = Z_tar(tar_ii)+(ts_ii-1)*dZ_tar(tar_ii);
        
%         % Signal Model in Polar System:
%         R_tar_ii = R_tar(tar_ii);
%         tau = 2*( R_tar_ii )/c;
%         cube_3d_ii = cube_3d_ii + exp( 1i*2*pi*( fc*tau + cr*t.*tau - cr*(tau.^2)/2 - ...
%                                                  2*X_rad2.*sind(theta_tar(tar_ii))/lambda - ...
%                                                  2*Z_rad2.*sind(ph_tar(tar_ii))/lambda ) );
        
        % Signal Model in Cartesian System:
        R_tar_ii = ...
            sqrt( ( X_tarii-X_rad2 ).^2 + ...
                  ( Y_tarii-Y_rad2 ).^2 + ...
                  ( Z_tarii-Z_rad2 ).^2 );

        tau = 2*( R_tar_ii )/c;   
        cube_3d_ii = cube_3d_ii + exp( 1i*2*pi*( fc*tau + cr*t.*tau - cr*(tau.^2)/2 ) );
    end
    cube_3dTS(:,:,:,ts_ii) = awgn( cube_3d_ii, snr );

end

%% ========================================================================
%% ========================================================================
%% ========================================================================
%% ========================================================================
