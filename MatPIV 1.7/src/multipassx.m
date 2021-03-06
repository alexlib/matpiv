function [x,y,u,v,SnR,Pkh]=multipassx(im1,im2,wins,Dt,overlap,sensit,maske,iter,datax,datay)

% MULTIPASSX - multiple passes
% function [x,y,u,v,snr,pkh]=multipassx(im1,im2,winsize,time,...
% overlap,sensit,maske,numofiterations,ustart,vstart)
%
% PIV in multiple passes to eliminate the displacement bias.
% Utilizes the increase in S/N by  halving the size of the
% interrogation windows after the first pass.
% Sub-function to MATPIV.
%
% See also:
%          MATPIV, SNRFILT, LOCALFILT, GLOBFILT, DEFINEWOCO

% Copyright 1998-2011, Kristian Sveen, jks@math.uio.no
% for use with MatPIV 1.7 and subsequent versions
% Distributed under the terms of the Gnu General Public License manager
% Time stamp: 09:20, Mar 4 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Image read
if ischar(im1)
    A=readmyimage(im1);
    B=readmyimage(im2);
else
    A=im1; B=im2;
end

A=double(A); B=double(B);

% background removal added 23. Oct 2002 - testversion
if exist([pwd,'/background.jpg'],'file')==2
    bgim=double(imread('background.jpg'));
    A=A-bgim; B=B-bgim;
    disp('   --> Background removed')
end
% Change October 1, 2002 - added simple image preconditioner:
%[A2,B2]=precondition(A,B);
% Check if preconditioning has made A,B==0. This may happen
% when the images are identical.
%if sum(A2(:))~=0 | sum(B2)~=0
%  A=A2; B=B2; clear A2 B2
%end

%%%%%%%%% First pass to estimate displacement in integer values:
if nargin==6
    maske=''; iter=3;
end
[sy,sx]=size(A);

if size(wins,1)==1
    if size(wins,2)==1
        wins=[wins, wins];
    end
    wins=repmat(wins,iter-2,1);
    wins=[wins; wins(end,:)/2];
    %note that wins is iter-1 long. The final iteration will be with the
    %same winsize as iter-1. In this case the windows are halved before
    %starting on the normalized correlation passes. It may be wiser to
    %complete another iteration with phase correlations first.
end

% check if datax and datay matrixes have correct size
if ~isempty(datax)
    if length(1:wins(1,1)*(1-overlap):sx-wins(1,1)+1)~=size(datax,2) || ...
            length(1:wins(1,2)*(1-overlap):sy-wins(1,2)+1)~=size(datax,1)
        disp(['  ===> Using supplied velocity field as input...with size ',num2str(size(datax)),'.'])
        [datax,datay]=velinterp(datax,datay,wins(1,:),overlap,size(A));
        disp(['Interpolated to ',num2str(size(datax))])
    else
        disp('No interpolation needed')
    end
end

for i=1:iter-1
    disp(['* Pass No: ',num2str(i)])
    if i==1 && ~isempty(datax)
        [x,y,datax,datay]=firstpass(A,B,wins(i,:),overlap,[],[],maske);
    else
        [x,y,datax,datay]=firstpass(A,B,wins(i,:),overlap,datax,datay,maske);
    end
    [datax,datay]=globfilt(x,y,datax,datay,3);
    [datax,datay]=localfilt(x,y,datax,datay,sensit,'median',3,maske);
    [datax,datay]=naninterp(datax,datay,'linear',maske,x,y);
    datax=floor(datax); datay=floor(datay);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % expand the velocity data to twice the original size
    if i~=iter-1
        if wins(i,1)~=wins(i+1,1)
            X=(1:((1-overlap)*2*wins(i+1,1)):sx-2*wins(i+1,1)+1) + wins(i+1,1);
            XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2;
        else
            XI=(1:((1-overlap)*wins(i+1,1)):sx-wins(i+1,1)+1)+(wins(i+1,1))/2; X=XI;
        end
        if wins(i,2)~=wins(i+1,2)
            Y=(1:((1-overlap)*2*wins(i+1,2)):sy-2*wins(i+1,2)+1) + wins(i+1,2);
            YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2;
        else
            YI=(1:((1-overlap)*wins(i+1,2)):sy-wins(i+1,2)+1)+(wins(i+1,2))/2; Y=YI;
        end
        disp('   Expanding velocity-field for next pass')
        datax=round(interp2(X,Y',datax,XI,YI'));
        datay=round(interp2(X,Y',datay,XI,YI'));
        [datax,datay]=naninterp(datax,datay,'linear',maske,...
            repmat(XI,size(datax,1),1),repmat(YI',1,size(datax,2)));
        datax=round(datax); datay=round(datay);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final pass. Gives displacement to subpixel accuracy.
disp('* Final Pass')
[x,y,u,v,SnR,Pkh]=finalpass(A,B,wins(end,:),overlap,round(datax),...
    round(datay),Dt,maske);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
