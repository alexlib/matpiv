% Example for PTV settings-file
% this file contains the necessery parameters to run MatPTV


% Basic parameters
ims={'mpim1b.bmp','mpim1c.bmp'}; % Images
%thr=[0.2 0.25 0.3 0.4 0.5 0.6]; % Threshold levels
thr=[0.6 0.9];
dt=0.04; %time between images


% Display settings
showfigs='yes'; % Wether to show particle positions during localization or not. 
                % It should save some time if you turn this off.
comlen=8; % length of comet trajectories to display during particle location



% particle property parameters
pmax=40; % maximum size of particles
pmin=3; % minimum size of particles





%particle matching parameters
maxdist=2; %Minimum distance between centroids of two particles
           %before they are considered to be the same.

	   
	   
	   
% Settings related to initial PIV guess

  % wether to use PIV or not to obtain a first guess at
  % the particle displacement. 'yes','no','onlyfirsttime'
  usepiv='onlyfirsttime'; 
  	      
  pivwin=64; % interrogation window size
  %dt is given above
  pivol=0.5; % overlap (in percent) of interrogation windows
  pivmet='mqd'; % calculation method. 'single','multi','multin','mqd'...
  pivwoc=''; % World coordinate file (consider deprecating this)
  pivmask=''; % Mask file
  %filter settings
  snrthr=1.3; %Signal to noise ration limit
  locthr=2.0; %Local filter limit
  locmet='median'; %Local filter method, 'median','mean'
  locker=3; %Size of the kernel in the local filter
  globthr=3; %Global filter threshold
  intermet='linear'; %Interpolation method, 'linear','weighted'