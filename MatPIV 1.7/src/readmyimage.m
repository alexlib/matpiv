function [out]=readmyimage(im)
% READMYIMAGE - helper function for reading images into MatPIV
%
% image=readmyimage('filename');
% This function is a very basic wrapper-function for reading images for use
% with MatPIV and MatPTV. The function basically removes the need for the
% Image Processing Toolbox
%
% Example: im=readmyimage('im04.bmp');
%
% SEE ALSO: MRGB2GRAY, RGB2GRAY, IND2GRAY, IND2RGB

imf=imfinfo(im);
[a,p]=imread(im);
wasuint8=isa(a,'uint8');
wasuint16=isa(a,'uint16');

if strcmp(imf.ColorType,'indexed')
    out=255*rgbconvert(ind2rgb(a,p));
elseif strcmp(imf.ColorType,'truecolor') %RGB images
    out=255*rgbconvert(a);
elseif strcmp(imf.ColorType,'grayscale')
    out=a;
else
    out=[];
    disp('Image color type not recognized');
end
if wasuint8 % change back to uint8 if that was the input
  out=uint8(out);
elseif wasuint16
  out=uint8(out);    
end

function [im]=rgbconvert(im)
% Helper code copied from mrgb2gray
% Returns a double, regardless of input
im=double(im);
[sx,sy,sz]=size(im);    
im=reshape(im,sx*sy,sz); %put image in an M*N by 3 array
im=im*[0.3 0.59 0.11]'; % multiply the weights, 0.3*R, 0.59*G, 0.11*B
im=reshape(im,sx,sy); % reshape the image back to correct size

