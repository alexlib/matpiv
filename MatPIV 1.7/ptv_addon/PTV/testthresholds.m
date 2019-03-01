function [b]=testthresholds(A,thr,act,centrdist,edgdist)
% TESTTHRESHOLDS - test a given list of threholds levels and count the
% total number of particles found.
%
% testthresholds(image,threholdslist);
%
% example: testthresholds('im04.bmp',[0.25:0.05:0.5]);
% the results can also be plotted:
% example: testthresholds('im04.bmp',[0.25:0.05:0.5],'on');
%
% minimum distance between newblob and oldblob centroids for 
% newblob to be accepted:
% testthresholds('im04.bmp',[0.25:0.05:0.5],'on',1);
%
% minimum distance between newblob and oldblob edges for 
% newblob to be accepted:
% testthresholds('im04.bmp',[0.25:0.05:0.5],'on',1,1);
%
if nargin<2, disp('Too few input arguments'); 
elseif nargin==2, act='on'; end 
if nargin==3, centrdist=1; edgdist=1; end
if nargin==4, edgdist=1; end

if ischar(A), 
  [A,p]=imread(A); 
  if misrgb(A), A=mrgb2gray(A); end, 
  %if ~isempty(p), A=ind2gray(A,p); end
end  
%A=double(A);



properties.pmax=40; % maximum size of particles
properties.pmin=3; % minimum size of particles
properties.xminsize=2; % minimum size in x direction, usually smart to set>1 to
           % avoid 1 pixel wide particles
properties.yminsize=2; 

thr=sort(thr);
loops=thr(end:-1:1);
tel=1;
for i=loops
    if tel==1
        [b,nb]=blobrec(A,i,properties);
    else    
        [b,nb]=blobrec(A,i,properties,b,centrdist,edgdist);
    end 
    tel=tel+1;
end 

disp([num2str(length(b.area)),' particles/blobs found'])

if strcmp(act,'on')
  %colors={'b','g','r','c','m','y'}; jalla=length(colors);
  %if length(colors)<length(loops)
  %  for i=length(colors)+1:length(loops)
  %    colors{i}=colors{i-jalla};
  %  end
  %end
  figure, imshow(A,'InitialMagnification','fit'), hold on
  plot(b.centr(:,1),b.centr(:,2),'b.')
  for i=1:length(b.area)
    plot([b.bound(i,1) b.bound(i,1) b.bound(i,1)+b.bound(i,3) b.bound(i,1)+...
	  b.bound(i,3) b.bound(i,1)],...
	 [b.bound(i,2) b.bound(i,2)+b.bound(i,4) b.bound(i,2)+b.bound(i,4)...
	  b.bound(i,2) b.bound(i,2)],'r-');
  end
end