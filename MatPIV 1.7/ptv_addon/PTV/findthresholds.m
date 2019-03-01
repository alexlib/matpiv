function [thres,blobs,st,ii]=findthresholds(A,uplo)
% FINDTHRESHOLDS - determine "optimal" thresholds for locating blobs in input image A
%
% findthresholds('mpim1b.bmp'); %will determine the optimal threshold
% levels for the image mpim1b.bmp. 
%
% UNDER DEVELOPEMENT - 
% o consider using the isodata-algorithm for threshold detection
% o get rid of for loop- takes too much time...
%
%

if ischar(A), 
  [A,p]=imread(A); 
  if isrgb(A), A=rgb2gray(A); end, 
  if ~isempty(p), A=ind2gray(A,p); end
end  

A=double(A);
uppr=max(A(:)); % upper value in image- no point checking above this
lowr=min(A(:)); % lower value in image- no point checking below this

tel=1;
if nargin==2,
    if length(uplo)==2, loops=uplo(2):-0.05:uplo(1); else, loops=uplo(end:-1:1); end 
else    
    loops=min([255 uppr]/255):-0.01:max([0 lowr]/255);
end 
%A=uint8(A); 
thres=[]; figure, hold on
for i=loops
    if tel==1
        [b,nb]=blobrec(A,i,3,40);
    else    
        [b,nb]=blobrec(A,i,3,40,b);
    end 
    if nb~=0, thres=[thres i]; end  
    tel=tel+1;
    %plot(i,nb,'*'); drawnow
end 
blobs=b;
[st,ii]=dbstack;
if ii==1
    disp([num2str(length(blobs.area)),' particles found, using thresholds [',...
            num2str(thres),']'])   
    imagesc(A), colormap(gray), axis ij
    plot(blobs.centr(:,1),blobs.centr(:,2),'bs')
end 