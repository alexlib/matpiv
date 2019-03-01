% Example for PTV settings-file
% this file contains the necessery parameters to run MatPTV


% Basic parameters
%ims={'im00.bmp','im04.bmp','im08.bmp'}; % Images
ims={'im0500.gif','im0501.gif','im0502.gif','im0503.gif'...
    'im0504.gif','im0505.gif','im0506.gif','im0507.gif',...
    'im0508.gif','im0509.gif'}; % Images
%thr=[0.2 0.25 0.3 0.4 0.5 0.6]; % Threshold levels
thr=[0.15 0.2 0.3 0.4 0.5 0.6 0.8];
dt=0.04; %time between images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display settings
showfigs='no'; % Wether to show particle positions during localization or not. 
                % It should save some time if you turn this off.
comlen=8; % length of comet trajectories to display during particle location

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% particle property parameters
props.pmax=40; % maximum size of particles
props.pmin=3;  % minimum size of particles
props.xminsize=2; % minimum size in x direction, usually smart to set>1 to
		       % avoid 1 pixel wide particles
props.yminsize=2; %as above	   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%particle parameters - Don't change these unless you know
%what you are doing
maxdist=3; %max distance between borders of two blobs before "new" blob accepted
edgedist=1;
%Matching parameters - same rules apply
m=2; %Power of error in matching distance (default=2)
maxerror=2; %max squared error in position estimate (default=2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings related to initial PIV guess

  % wether to use PIV or not to obtain a first guess at
  % the particle displacement. 'yes','onlyfirsttime'
  usepiv='yes'; % default = 'onlyfirsttime'
  whodoyoutrust=[1 0]; % which to trust from the PTV and PIV
		       % measurement. This should be a 1x2 vector with
		       % sum==1. The number represent weights for PTV
		       % and PIV respectively. This means that a
		       % vector = [0.5 0.5] will use the average of
		       % the PIV and the PTV predictions, while
		       % whodoyoutrust=[1 0] means you only use the
		       % PTV predictions to estimate position in the
		       % next frame.
  pivwin=64; % interrogation window size
  %dt is given above
  pivol=0.5; % overlap (in percent) of interrogation windows
  pivmet='mqd'; % calculation method. 'single','multi','multin','mqd'...
  pivwoc=''; % World coordinate file (consider deprecating this)
  pivmask=''; % Mask file
  %filter settings
  snrthr=1.3; %Signal to noise ration limit
  locthr=2.5; %Local filter limit
  locmet='median'; %Local filter method, 'median','mean'
  locker=3; %Size of the kernel in the local filter
  globthr=3; %Global filter threshold
  intermet='linear'; %Interpolation method, 'linear','weighted'
  
  %{
  res = matptv('myptvsettings');
  figure
  for i = 2:length(res)
    x = res(i).blobs.centr(:,1);
    y = res(i).blobs.centr(:,2);
    u = res(i).blobs.ptvvel(:,1);
    v = res(i).blobs.ptvvel(:,2);
    plot(x,y,'.')
    hold on
    quiver(x,y,u,v,1)
    drawnow
    pause(.1)
    hold off
   end
  %}   
  
  