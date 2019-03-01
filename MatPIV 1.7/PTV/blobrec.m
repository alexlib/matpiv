function [blobs,numofnew]=blobrec(a,thr,lowlim,uplim,oldblobs,centrdist,edgdist)
% BLOBREC - Find blobs in a grayscale image
% function [blobs,numofnew]=blobrec(image,threshold,minsize,maxsize,oldblobs,maxdist)
%
% This function will recognise blobs in an input image a, with threshold 
% between 0 and 1. Optional input is the lower and upper limits on the 
% size/area (in pixels) of the recognised blobs (area; defaults are 3 and 40).
%
% If OLDBLOBS is included, all new blobs will be added to the
% OLDBLOBS-structure and returned.
%
% CENTRDIST
%
% EDGEDIST
%
% Copyright Jul 23 2000 - March 11 2011, jks@math.uio.no/matpiv@gmail.com

if nargin<2
    disp('Wrong input argument!'); return
elseif nargin<3
    lowlim=3; uplim=40;
end

if ischar(a), 
    a=readmyimage(a); 
%  if isrgb(a), a=rgb2gray(a); end
%  if ~isempty(p), a=ind2gray(a,p); end
  %a=double(a); 
end

if isstruct(lowlim)
  strustat=1;
  if nargin>3
    edgdist=centrdist;  centrdist=oldblobs;  oldblobs=uplim;  
  end
  uplim=lowlim.pmax; xmin=lowlim.xminsize;  
  ymin=lowlim.yminsize; lowlim=lowlim.pmin;
else
  xmin=2; ymin=2; strustat=0;
end
  
[sta,ind]=dbstack;

if nargin>=4
  xtrastring='new ';
else
  xtrastring='';
end
%b=im2bw(a,thr);     
b=(a > 255*thr); % this is essentially what im2bw does
l = bwlabel(b,4);
%b(b==0)=nan;
%l=label(b,4);

%regionprops is supposed to replace imfeature, but seems 25 times slower!
%(on Jan30 2004)
stats2 = regionprops(l,'basic'); %toc
%stats2 = imfeature(l,'basic');

area=cat(1,stats2.Area);
centr=cat(1,stats2.Centroid);
bounding=cat(1,stats2.BoundingBox);
tmpbx=bounding(:,3);
tmpby=bounding(:,4);

%Remove particles that are larger or smaller than our min/max sizes
%and additionally remove particles that are smaller then 2 pixels in
%either dimension.
centr=centr(area>=lowlim & area <= uplim &...
	    tmpbx>=xmin & tmpby>=ymin,:);
bounding=bounding(area>=lowlim & area <= uplim &...
		  tmpbx>=xmin & tmpby>=ymin,:); 
area=area(area>=lowlim & area <= uplim &...
	  tmpbx>=xmin & tmpby>=ymin); 


if nargin>4
  blobs.area=area; blobs.centr=centr; blobs.bound=bounding;
  blobs=compareblobs(oldblobs,blobs,centrdist,edgdist);
  numofnew=length(blobs.area)-length(oldblobs.area);
else    
  blobs.area=area; blobs.centr=centr; blobs.bound=bounding;
  numofnew=length(area);
end 

% Output of result to screen...
% We do not want this message printed at each threhold level for all images
% in the case where we run this from MATPTV. Therefore suppress output...
if size(sta,1)<=1
    disp(['Number of ',xtrastring,'blobs: ', num2str(numofnew), ...
            '. Threshold used ', num2str(thr)]);
else
    if isempty(findstr(sta(end).name,'matptv'))   
        disp(['Number of ',xtrastring,'blobs: ', num2str(numofnew), ...
                '. Threshold used ', num2str(thr)]);
    end
end 
