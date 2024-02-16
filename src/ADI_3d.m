function adi = ADI_3d(SLC_TS3d)
% SLC_TS: r*az*el*t
slc_ts_std = std( abs(SLC_TS3d),0,4 );
slc_ts_mean = mean( abs(SLC_TS3d),4 );
adi = slc_ts_std./slc_ts_mean;
end