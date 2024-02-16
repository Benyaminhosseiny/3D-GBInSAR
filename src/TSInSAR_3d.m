function [TS_phase,TS_intf_phase,TS_cum_phase] = TSInSAR_3d(TS_SLC_3d)
% TS_SLC: r_dir x az_dir x el_dir x t

TS_phase      = angle( TS_SLC_3d );
TS_intf_phase = wrapToPi( diff(TS_phase,1,4) );
TS_cum_phase  = cumsum(TS_intf_phase,4);

end