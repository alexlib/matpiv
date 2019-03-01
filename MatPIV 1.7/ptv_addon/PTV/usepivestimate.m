function [pivout]=usepivestimate(ims1,ims2,setfile)
%
%
%
%
eval(setfile)

% 
pivout=matpiv(ims1,ims2,pivwin,dt,pivol,pivmet,pivwoc,pivmask);

%filtering
[su,sv]=snrfilt(pivout.x,pivout.y,pivout.u,pivout.v,pivout.snr,snrthr);

[mu,mv]=localfilt(pivout.x,pivout.y,su,sv,locthr,locmet,locker,pivmask);

[gu,gv]=globfilt(pivout.x,pivout.y,mu,mv,globthr);

[fu,fv]=naninterp(gu,gv,intermet,pivmask,pivout.x,pivout.y);

%allocate the filtered velocities to the output structure
pivout.fv=fv; pivout.fu=fu;