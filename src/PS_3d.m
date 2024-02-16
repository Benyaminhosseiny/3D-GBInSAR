function ps_mask = PS_3d(SLC_TS_3d,thresh_adi,thresh_mag)
% PS detection:
% SLC_TS: r*c*t
% thresh_mag from max: dB

if nargin<3
    thresh_mag = 15; %dB
end
if nargin<2
    thresh_adi = 0.1;
end

adi = ADI_3d(SLC_TS_3d);
slc_abs = mean( 10*log10(abs(SLC_TS_3d)),4 ); % dB
% slc_abs = 10*log10(abs(SLC_TS_3d(:,:,:,end))); % dB
%
thresh_mag = max(slc_abs,[],'all')-thresh_mag; %dB
mask_adi = adi; mask_adi(adi>thresh_adi) = 0; mask_adi(adi<=thresh_adi) = 1;
mask_mag = slc_abs; mask_mag(slc_abs<thresh_mag) = 0;mask_mag(slc_abs>=thresh_mag) = 1;
ps_mask = mask_adi.*mask_mag;

end