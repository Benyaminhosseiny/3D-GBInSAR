function [d_hat_3D_cartesian_dif,d_hat_3D_cartesian_ts,d_hat_3D_cartesian_total] = Displacement_vec_Cartesian(Tar3dLoc,SAR3dLoc,RC_sig_PS_ts,lambda,clutter_rmv_flag)
% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. "Structural displacement 
% monitoring using ground-based synthetic aperture radar." International Journal of 
% Applied Earth Observation and Geoinformation (JAG) 116 (2023): 103144.
% https://doi.org/10.1016/j.jag.2022.103144

% Tar3dLoc    : Target's 3D position in local coordinate system.            [3x1 vector]
% SAR3dLoc    : Radar's 3D position during data aquisition for SAR imaging. [3d array NxMx3 ]
% RC_sig_PS_ts: Range-compressed signals corresponding to the detected targets: 3d array [timeseries x targets x sar_steps]
% lambda      : Wavelength (m)
% clutter_rmv_flag: clutter removal flag 0 or 1.

% OUTPUTS:
% d_hat_3D_cartesian_dif  : amount of targets' differential displacemet vector at each time (epoch). 3d array: [epochs x targets x xyz]
% d_hat_3D_cartesian_ts   : cumulative displacemet vector at each time.                              3d array: [epochs x targets x xyz]
% d_hat_3D_cartesian_total: total displacemet vector after time-series.                              2d array: [targets x xyz]

num_PS      = size(Tar3dLoc,1);
radar_steps = size(RC_sig_PS_ts,3);
% Clutter removal:
if clutter_rmv_flag
    for ps_ii = 1:num_PS
        for radar_ii = 1:radar_steps
            circlefit_Par = CircleFitByPratt( [real(RC_sig_PS_ts(:,ps_ii,radar_ii)),imag(RC_sig_PS_ts(:,ps_ii,radar_ii))] );
            x_c = circlefit_Par(1);
            y_c = circlefit_Par(2);
            radius_c = circlefit_Par(3);
            RC_sig_PS_ts(:,ps_ii,radar_ii) = RC_sig_PS_ts(:,ps_ii,radar_ii)-(x_c+1i*y_c);
        end
    end
end
% Phase only:
RC_angle_sig_PS_ts=angle(RC_sig_PS_ts); % 3d array: [time_samples x num_PS x radar_steps]



defo_rc_PS = wrapToPi( diff(RC_angle_sig_PS_ts,1) ) * lambda / 4 / pi; % 3d array [epochs x targets x sar steps]

epochs = size(defo_rc_PS,1);
for epoch_i = 1:epochs
    for ps_i = 1:num_PS
        Tar3dLocii = Tar3dLoc(ps_i,:);
        
        d_hat_3D_cartesian_dif(epoch_i,ps_i,:) = Dhat3D_vec_Cartesian( Tar3dLocii,SAR3dLoc,squeeze(defo_rc_PS(epoch_i,ps_i,:)),lambda );
        
        % New target locations after 3D displacement estimation:
        Tar3dLoc(ps_i,:) = Tar3dLoc(ps_i,:)+squeeze(d_hat_3D_cartesian_dif(epoch_i,ps_i,:))';
    end
    
end
d_hat_3D_cartesian_ts = cumsum(d_hat_3D_cartesian_dif,1);

d_hat_3D_cartesian_total = squeeze(d_hat_3D_cartesian_ts(end,:,:)); % targets x xyz
if num_PS==1
    d_hat_3D_cartesian_total = d_hat_3D_cartesian_total';
end



%%
    function d_hat = Dhat3D_vec_Cartesian(Tar3dLoc,SAR3dLoc,defo_rc_target,lambda)
        % Tar3dLoc: Target's 3D position in local coordinate system. [3x1 vector]
        % SAR3dLoc: Radar's 3D position during data aquisition for SAR imaging. [NxMx3 array]
        % defo_rc_target: all the LOS interferometric measurements of range-compressed signal corresponding to the target [1d vector]
        
        % A 2020 3-D Deformation Measurement Based on Three GB-MIMO Radar Systems: Experimental Verification and Accuracy Analysis
        
        Tar3dLoc = Tar3dLoc(:);
        
        x_n = SAR3dLoc(:,:,1);
        y_n = SAR3dLoc(:,:,2);
        z_n = SAR3dLoc(:,:,3);
        radar_position = [x_n(:), y_n(:), z_n(:)]';
        defo_rc_target = defo_rc_target(:);
        
        
        U_L = Tar3dLoc-radar_position;
        U_L = U_L';
        U_L = U_L./( sqrt( U_L(:,1).^2+U_L(:,2).^2+U_L(:,3).^2 ) );
        d_hat=lsqlin( U_L, defo_rc_target );
        % d_hat=lsqlin( U_L, defo_rc_target, ones(size(U_L)),(lambda/4)*ones(size(defo_rc_target)) ); % inequality costraint 1*x<lambda/4
%         d_hat=lsqlin( U_L, defo_rc_target,[],[],[],[],-(lambda/4)*ones(3,1),(lambda/4)*ones(3,1) ); % upper and lower bound costraint 
%         d_hat=lsqlin( U_L, defo_rc_target,[],[],[],[],[0,0,-lambda/4],[0,0,lambda/4] );

        % Geometry costraint:
        theta = atan( Tar3dLoc(2)/Tar3dLoc(1) );
        phi   = acot( Tar3dLoc(2)/Tar3dLoc(3)/sin(theta) );
        T     = [sin(theta);sin(phi);sqrt(1-sin(theta).^2-sin(phi).^2)].^(-1);
        T(1)=0;
%         d_hat=lsqlin( U_L, defo_rc_target,[],[],[],[],-(lambda/4)*abs(T),(lambda/4)*abs(T) ); % upper and lower bound costraint
        
    end

end