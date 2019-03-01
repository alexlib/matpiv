function [x,y,u,v]=pixel2world(upix,vpix,xpix,ypix,lswo1,lswo2,mapping)
% PIXEL2WORLD - transform from pixels and pixels/sec to centimeters and cm/s
%
% [x,y,u,v]=pixel2world(upix,vpix,xpix,ypix,lswo1,lswo2,mapping);
% alternatively:
% [x,y,u,v]=pixel2world(u,v,x,y,'worldco1.mat','linear');
%
% Calculates the pixels to world coordinate transformation
% You need first to use the m-file DEFINEWOCO to specify your
% world coordinate system. Definewoco then calculates 6 numbers
% that are saved to file 
%
% UPIX, VPIX, XPIX, YPIX are the velocity field in pixels and pixels/second
% 
% LSWO1 and LSWO2 are the mapping factors for x and y direction
% respectively.
%
% MAPPING is the mapping function from pixel to world coordinates
% and can be any of 
%                         'linear'    - fits the function 1 + x + y
%                         'nonlinear' - fits 1+x+y+x^2+y^2+xy
%                         'bilinear'  - fits 1+x+y+xy
% Notice that the actual function fitting is performed inside DEFINEWOCO so
% you need to apply the same MAPPING here as you did when you defined the
% coordinate system.
%
% Contributions from Niclas Hjerdt
%
% See also: DEFINEWOCO, MATPIV

% 1998 - 2005 , jks@math.uio.no
% For use with MatPIV 1.6 and subsequent versions
%
% Copyright J.K.Sveen (jks@math.uio.no)
% Dept. of Mathematics, Mechanics Division, University of Oslo, Norway
% Distributed under the terms of the Gnu General Public License

if nargin==6 
    mapping='linear';  
    if ischar(lswo1)
        l=load(lswo1);
        lswo1=l.comap(:,1);
        lswo2=l.comap(:,2);
    end
elseif nargin==5
    mapping = 'linear';
    if ischar(lswo1)
        l=load(lswo1);
        lswo1=l.comap(:,1);
        lswo2=l.comap(:,2);
    else
        lswo2=lswo1(:,1);
        lswo1=lswo1(:,2);   
    end
end


if strcmp(mapping,'linear')==1
  lswo1(4:6)=0; lswo2(4:6)=0;
elseif strcmp(mapping,'nonlinear')==1 | strcmp(mapping,'bilinear')==1
  if length(lswo1)<4 
    disp('This mapping file only contains enough defined points for a linear mapping.');
    disp('You need to redefine your world coordinate points using DEFINEWOCO.M')
    return
elseif length(lswo1)==4
    lswo1(5:6)=0; lswo2(5:6)=0;    
end
elseif strcmp(mapping,'nonlinear')==0 & strcmp(mapping,'linear')==0 & strcmp(mapping,'bilinear')==0
  disp('No such mapping available, try `linear` or `nonlinear`!')
  return
end

fprintf(['* Calculating the pixel to world transformation using ',mapping,' mapping'])
for ii=1:1:size(upix,2)
  for jj=1:1:size(upix,1)
    u(jj,ii)=(lswo1(2)*upix(jj,ii)) + (lswo1(3)*vpix(jj,ii));
    v(jj,ii)=(lswo2(2)*upix(jj,ii)) + (lswo2(3)*vpix(jj,ii));
    x(jj,ii)=lswo1(1)+ lswo1(2)*xpix(jj,ii)+ lswo1(3)*ypix(jj,ii)+...
	  lswo1(4)*(xpix(jj,ii).*ypix(jj,ii))+...
	  lswo1(5)*(xpix(jj,ii).^2)+lswo1(6)*(ypix(jj,ii).^2);
    y(jj,ii)=lswo2(1)+ lswo2(2)*xpix(jj,ii)+ lswo2(3)*ypix(jj,ii)+...
	  lswo2(4)*(xpix(jj,ii).*ypix(jj,ii))+...
	  lswo2(5)*(xpix(jj,ii).^2)+lswo2(6)*(ypix(jj,ii).^2);
  end
end

fprintf(' - DONE\n')
