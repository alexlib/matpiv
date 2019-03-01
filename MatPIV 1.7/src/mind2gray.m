function [fname]=mind2gray(fname,cmap);
% MIND2GRAY - convert indexed image to grayscale
%
% This is a cheap hack intended for users with no Image Processing Toolbox
% license. 
%
%

if nargin==1 & ischar(fname)
    [fname,cmap]=imread(fname);
end
[sx,sy]=size(fname);
cmaplength=size(cmap,1);
% fname=double(fname(:));
% cmap=sum(cmap')';%cmap(:,1);
% 
% for i=1:cmaplength
% %    ind=find(fname==i);
%     fname(fname==i)=fname(fname==i)*cmap(i);
% end
% 
% fname=reshape(fname,sx,sy);

imwrite(fname,cmap,'test.jpg','JPEG','Quality',100); 
[fname,emptycmap]=imread('test.jpg'); 
delete test.jpg
fname=mrgb2gray(fname);
%figure
%h=image(fname);
%colormap(cmap)
%h=get(gca);
%h=get(h.Children);
%fname=h.CData;
%close