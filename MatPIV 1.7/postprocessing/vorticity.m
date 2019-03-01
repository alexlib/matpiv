function [omega,str]=vorticity(x,y,u,v,method)

% [omega,str]=vorticity(x,y,u,v,method)
% 
% This function calculates differential quantities from a given
% flowfield  x y u v. METHOD should be one of 'circulation', 
% 'richardson', 'leastsq' or 'centered'. Default is least squares.
% 'centered' uses the MATLAB function CURL.
%
% Notice that the 'richardson' and 'leastsq' options are based on a
% regular axis coordinate system (i.e. NOT image axis which is
% reversed in the vertical direction). This means that for
% measurements in an image frame of reference (i.e. NOT converted to
% world coordinates) you may have to input -v instead of v.
%
% The second output is str=0.5*(du/dy+dv/dx)
%
% EXAMPLES: 
%         w=vorticity(x,y,u,v); pcolor(x(3:end-2,3:end-2),y(3:end-2,3:end-2),w);
%         w=vorticity(x,y,u,v,'circulation'); ...
%           pcolor(x(2:end-1,2:end-1),y(2:end-1,2:end-1),w);
%
%   To test the accuracy of MatPIV and vorticity:
%      [im1,im2,U]=makeimage3d(1024,1024,30000,'rankine');
%      Wi=interp2(U.x,U.y,U.vort,x,y); % synthetic vorticity field
%      [x,y,u,v,snr]=matpiv(im1,im2,[32 32;32 32],1,0.5,'multinfft');
%      [su,sv]=snrfilt(x,y,u,v,snr,2);
%      [fu,fv]=naninterp(su,sv);
%      w1=vorticity(x,y,fu,fv,'leastsq');
%      w2=vorticity(x,y,fu,fv,'circulation');
%      w3=vorticity(x,y,fu,fv,'centered');
%      w4=vorticity(x,y,fu,fv,'richardson');
%     figure, 
%      subplot(2,2,1), pcolor(x(3:end-2,3:end-2),y(3:end-2,3:end-2),w1./max(Wi(:)))
%      caxis([0 1]), title('Least Squares')
%      subplot(2,2,2),
%      pcolor(x(2:end-1,2:end-1),y(2:end-1,2:end-1),w2./max(Wi(:)))
%      caxis([0 1]), title('Circulation')
%      subplot(2,2,3), pcolor(x,y,w3./max(Wi(:)))
%      caxis([0 1]), title('Centered (curl.m)')
%      subplot(2,2,4),pcolor(x(3:end-2,3:end-2),y(3:end-2,3:end-2),w4./max(Wi(:)))
%      caxis([0 1]), title('Richardson')
%
%   Notice that the vorticity in this example has been normalized by the
%   max vorticity in the synthetic velocity field so that the
%   color-scale is between 0 and 1.
%
%
%  See also: MATPIV, MAKEIMAGE3D, VORTICITY, CURL, INTERP2, 


% Oct 11, 2005
% Copyright Kristian Sveen, jks@math.uio.no
% for use with MatPIV 1.6.+
if nargin ==4 | nargin==1
  method='leastsq';
end
if nargin==2 | nargin==1
  if ischar(x)
    if nargin==1
      method='leastsq';
    else
      method=y;
    end
    vel=load(x);
    x=vel(:,1); y=vel(:,2); u=vel(:,3); v=vel(:,4);
  else
    disp('Wrong input to VORTICITY');
    return
  end
end

if size(x,2)==1
  disp('Converting vectors to matrices')
  [x,y,u,v]=fixdigim(x,y,u,v);
end




% SCALE is the scale for velocity vectors.
scale=2/max(sqrt(u(:).^2 + v(:).^2));

DeltaX=x(1,2)-x(1,1);
DeltaY=y(1,1)-y(2,1);

if strcmp(method,'circulation')==1
  Vort=zeros(size(x)-2);
  str=Vort;
  Dx=(1/(8*DeltaX)).*[-1 0 1
    -2 0 2
    -1 0 1];
  Dy=(1/(8*DeltaY)).*[1 2 1
    0 0 0
    -1 -2 -1];  
  Vort=conv2(v,Dx,'valid')-conv2(u,Dy,'valid');
  outp=-real(Vort);
  str=0.5*(conv2(v,Dx,'valid')+conv2(u,Dy,'valid'));
  xa=x(2:end-1,2:end-1);
  ya=y(2:end-1,2:end-1);
  
elseif strcmp(method,'richardson')==1
  vor=zeros(size(x)-4);
  str=vor;
  for i=3:1:size(x,2)-2
    for j=3:1:size(x,1)-2
      tmp1= -(-v(j,i+2) +8*v(j,i+1) -8*v(j,i-1) +...
	      v(j,i-2))/(12*DeltaX);
      tmp2=+ (-u(j+2,i) +8*u(j+1,i) -8*u(j-1,i) +...
	      u(j-2,i))/(12*DeltaY);
      vor(j-2,i-2)= tmp1+tmp2;
      str(j-2,i-2)=0.5*(-tmp1+tmp2);
    end
  end
  outp=vor;
  xa=x(3:end-2,3:end-2);
  ya=y(3:end-2,3:end-2);
  
elseif strcmp(method,'leastsq')==1
  vor=zeros(size(x)-4);
  str=vor;
  for i=3:1:size(x,2)-2
    for j=3:1:size(x,1)-2
      tmp1=-(2*v(j,i+2) +v(j,i+1) -v(j,i-1) -2*v(j,i-2))/(10*DeltaX);
      tmp2=+(2*u(j+2,i) +u(j+1,i) -u(j-1,i) -2*u(j-2,i))/(10*DeltaY);
      vor(j-2,i-2)= tmp1+tmp2;
      str(j-2,i-2)= 0.5*(-tmp1+tmp2);
	    
    end
  end
  outp=vor;
  xa=x(3:end-2,3:end-2);
  ya=y(3:end-2,3:end-2);
elseif strcmp(method,'centered')==1 & exist('curl','file')==2
  [vor,str]=curl(x,y,u,v); % str=angular momentum here
  outp=vor;
  xa=x(1,:);ya=y(:,1);
else
  disp([method,' is not a valid calculation method!!!'])
  disp(['Check your spelling or that you have the function ',...
	'CURL in your path']);
  return
end

if nargout==0
    %plot the stuff if no output argument is specified
    figure
    pcolor(xa,ya,outp);
    hold on
    shading interp
    vekplot2(x(:).',y(:).',u(:).',v(:).',scale,'k');
    axis tight
elseif nargout>=1
    omega=outp;
end
