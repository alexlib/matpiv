%% simplest

[x,y,u,v]=matpiv('im00.bmp','im04.bmp',64,0.0012,0.5,'single');
quiver(x,y,u,v)

%% extensive

% definewoco('mpwoco.bmp','.');
definewoco('woco.bmp','o');
mask('mpim1b.bmp','worldco.mat');
[x,y,u,v,snr,pkh]=matpiv('im00.bmp','im04.bmp',...
[64 64;64 64;32 32;32 32;16 16;16 16],...
0.0012,0.5,'multin','worldco.mat','polymask.mat');
[su,sv]=snrfilt(x,y,u,v,snr,1.3);
[pu,pv]=peakfilt(x,y,su,sv,pkh,0.5);
[gu,gv]=globfilt(x,y,pu,pv,3);
[mu,mv]=localfilt(x,y,gu,gv,2,'median',3,'polymask.mat');
[fu,fv]=naninterp(mu,mv,'linear','polymask.mat',x,y);

% results
quiver(x(1:4:end),y(1:4:end),fu(1:4:end),fv(1:4:end),2); 
axis tight

w=magnitude(x,y,fu,fv);
pcolor(x,y,w)
shading flat
colorbar
w=vorticity(x,y,fu,fv,'circulation');
pcolor(x(2:end-1,2:end-1),y(2:end-1,2:end-1),w)
shading flat
colorbar
h=mstreamline(x,y,fu,fv,2);