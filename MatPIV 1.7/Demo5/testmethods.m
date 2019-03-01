
im1='cyl2.jpg';
im2='cyl3.jpg';
win=64;
ol=0.75;
sthr=1.2;
gthr=1;
lthr=2.5;

[x1,y1,u1,v1,snr1,pkh1]=matpiv(im1,im2,win,0.001,ol,...
'phase','worldco2.mat','polymask2.mat');
[x2,y2,u2,v2,snr2,pkh2]=matpiv(im1,im2,win,0.001,ol,...
'single','worldco2.mat','polymask2.mat');
[x3,y3,u3,v3,snr3,pkh3]=matpiv(im1,im2,win,0.001,ol,...
'mqd','worldco2.mat','polymask2.mat');
[x4,y4,u4,v4,snr4,pkh4]=matpiv(im1,im2,win,0.001,ol,...
'norm','worldco2.mat','polymask2.mat');

%scale=3/nanmax(sqrt(u1(:).^2 + v1(:).^2));


[su,sv]=snrfilt(x1,y1,u1,v1,snr1,sthr);
[su,sv]=globfilt(x1,y1,su,sv,gthr);
[su1,sv1]=localfilt(x1,y1,su,sv,lthr,'median');
[fu1,fv1]=naninterp(su1,sv1,'linear','polymask2.mat',x1,y1);

[su,sv]=snrfilt(x2,y2,u2,v2,snr2,sthr);
[su,sv]=globfilt(x2,y2,su,sv,gthr);
[su2,sv2]=localfilt(x2,y2,su,sv,lthr,'median');
[fu2,fv2]=naninterp(su2,sv2,'linear','polymask2.mat',x2,y2);

[su,sv]=snrfilt(x3,y3,u3,v3,snr3,sthr);
[su,sv]=globfilt(x3,y3,su,sv,gthr);
[su3,sv3]=localfilt(x3,y3,su,sv,lthr,'median');
[fu3,fv3]=naninterp(su3,sv3,'linear','polymask2.mat',x3,y3);

[su,sv]=snrfilt(x4,y4,u4,v4,snr4,sthr);
[su,sv]=globfilt(x4,y4,su,sv,gthr);
[su4,sv4]=localfilt(x4,y4,su,sv,lthr,'median');
[fu4,fv4]=naninterp(su4,sv4,'linear','polymask2.mat',x4,y4);


scale1=1/nanmax(sqrt(u1(:).^2 + v1(:).^2));
scale2=1/nanmax(sqrt(su1(:).^2 + sv1(:).^2));
scale3=1/nanmax(sqrt(fu1(:).^2 + fv1(:).^2));

figure
subplot(2,2,1)
vekplot2(x1,y1,u1,v1,scale1), axis([-5 20 28 54])
subplot(2,2,2)
vekplot2(x2,y2,u2,v2,scale1), axis([-5 20 28 54])
subplot(2,2,3)
vekplot2(x3,y3,u3,v3,scale1), axis([-5 20 28 54])
subplot(2,2,4)
vekplot2(x4,y4,u4,v4,scale1), axis([-5 20 28 54])
%subtitle('Unfiltered')


figure
subplot(2,2,1)
vekplot2(x1,y1,su1,sv1,scale2), axis([-5 20 28 54])
subplot(2,2,2)
vekplot2(x2,y2,su2,sv2,scale2), axis([-5 20 28 54])
subplot(2,2,3)
vekplot2(x3,y3,su3,sv3,scale2), axis([-5 20 28 54])
subplot(2,2,4)
vekplot2(x4,y4,su4,sv4,scale2), axis([-5 20 28 54])
%subtitle('Filtered')

figure
subplot(2,2,1)
vekplot2(x1,y1,fu1,fv1,scale3), axis([-5 20 28 54])
subplot(2,2,2)
vekplot2(x2,y2,fu2,fv2,scale3), axis([-5 20 28 54])
subplot(2,2,3)
vekplot2(x3,y3,fu3,fv3,scale3), axis([-5 20 28 54])
subplot(2,2,4)
vekplot2(x4,y4,fu4,fv4,scale3), axis([-5 20 28 54])
%subtitle('Filtered and interpolated')
