function [out]=misrgb(im);
% MISRGB - checks if input is an RGB image or not
%
% if the input is of type 'double', the pixel values need to be 0<=P<=1
%
%
% output is logical

% v 0.2
% J.K.Sveen@damtp.cam.ac.uk, 2004, August 27.
% Distributed under the terms of the GNU GPL:
% http://www.gnu.org/copyleft/gpl.html

if ischar(im)
  im=imread(im);
end

[sx,sy,sz]=size(im);
out=0;


if sz==3
  out=1;
  if isa(im,'double') % if the image is a double, we need to check
                      % its values
    mmax=max(im(:));
    mmin=min(im(:));
    if ~(mmin>=0 & mmax<=1) % if the values are outside, the image is
                            % not considered to be RGB
      out=0;
    end
  end
end

out=logical(out);