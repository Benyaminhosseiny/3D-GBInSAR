function [rc,slc3d] = create_slc_3d(cube_3d,r_fft,az_fft,el_fft)
% cube_3d: range_dir x az_dir x el_dir

rc  = fft( cube_3d,r_fft,1 );
azc = fftshift( fft( rc,az_fft,2 ),2 );
slc3d = fftshift( fft( azc,el_fft,3 ),3 );
end