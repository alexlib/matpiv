function [u,v,info1]=unifilter(x,y,u,v,thr,msk);
% UNIFILTER - Universal outlier detection for PIV data
%
% Shamelessly copied from Westerweel and Scarano, Exp Fluids, 2005
%
% [u,v]=unifilter(x,y,u,v,thr,msk);
% 
% x,y,u,v - velocity arrays
% thr - normalized threshold (default 2)
% msk - MatPIV mask file (default none)
%
% example: 
%
%        [x,y,u,v,snr]=matpiv('mpim1b.bmp','mpim1c.bmp',[64 64;32 32;32 32],...
%                      0.008,0.5,'multinfft');
%        [su,sv]=unifilt(x,y,u,v);
%
% NEED TO IMPLEMENT MASKING

if nargin==2
  thr=2; u=x; v=y; msk='';
elseif nargin==3
  u=x; v=y; msk='';
elseif nargin==4
  if ischar(v)
    thr=u; msk=v;
    u=x; v=y; 
  else
    thr=2; msk='';
  end
elseif nargin==5
  msk='';
end

if ~isempty(msk)
  msk=load(msk);
end

[sy,sx]=size(u);

Medianres=zeros(sy,sx);
NormFluct = zeros(sy,sx,2);
b=1;
eps=0.1;

for c=1:2
  if c==1; VelComp=u; else, VelComp=v; end
  for i=1+b:sx-b
    for j=1+b:sy-b
      Neigh=VelComp(j-b:j+b,i-b:i+b);
      NeighCol=Neigh(:);
      NeighCol2=[NeighCol(1:(2*b+1)*b+b); NeighCol((2*b+1)*b+b+2:end)];
      Median=mnanmedian(NeighCol2);
      Fluct=VelComp(j,i)-Median;
      Res=NeighCol2-Median;
      MedianRes=mnanmedian(abs(Res));
      NormFluct(j,i,c)=abs(Fluct/(MedianRes+eps));
    end
  end
end

info1=(sqrt(NormFluct(:,:,1).^2 + NormFluct(:,:,2).^2)>thr);
u(info1)=nan;
v(info1)=nan;
