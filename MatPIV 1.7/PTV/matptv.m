function [particles]=matptv(setfile)
% MATPTV - Particle Tracking Velocimetry addon to MatPIV
%
% particles=matptv(settingsfile);
% will track all particles in the IMAGES specified in SETTINGSFILE (either
% an .avi movie, or a comma separated list of images in a cell array).
% To do this MATPTV also uses a number of parameters specified in
% this file. SETTINGSFILE is an m-file (not a function) that can be
% executed within MATPTV to yield the desired settings.
%
%

% time stamp: 16:00, April 15 2014
% c, jks@math.uio.no
%
if nargin<1, disp('please specify settings-file'); return,
else
    eval(setfile)
    if ischar(ims)
        if strcmp(ims(end-3:end),'.avi'),
            ims=aviread(ims,aviframeno);
        else, disp('Only AVI movies are supported'); return, end
    elseif iscell(ims)
        for i=1:size(ims,2)
            ims2(i).cdata=readmyimage(strvcat(ims(:,i)));
            %[tmp,cm]=imread(strvcat(ims(:,i)));
            %if isrgb(tmp), tmp=rgb2gray(tmp); end
            %if ~isempty(cm), tmp=ind2gray(tmp,cm); end
            %ims2(i).cdata=tmp; ims2(i).colormap=cm;
        end
        ims=ims2; clear ims2
    end
end
thr=sort(thr);

comtel=1;
for i=1:size(ims,2)
    %run through each image
    disp([' - Locating particles in frame number ',num2str(i)])
    %need to implement a masking strategy - particles need to be born
    %and die at mask border - not at image edges
    %  if ~isempty(pivmask), msk=load(pivmask); end
    
    for th=length(thr):-1:1
        %now run through all threshold levels for this image
        if th==length(thr)
            clear blobs
            blobs=blobrec(ims(i).cdata,thr(th),pproperties);
        else
            blobs=blobrec(ims(i).cdata,thr(th),pproperties,blobs,maxdist,edgedist);
        end
        particles(i).t=(i-1)*dt;
        particles(i).blobs=blobs;
    end
    disp(['  -> Located ',num2str(length(blobs.area)),' particles'])
    
    if (strcmp(usepiv,'yes') & i~=size(ims,2)) | (strcmp(usepiv,'onlyfirsttime') & i==1)
        %here we compute an estimated velocity field using PIV. This
        %field will subsequently be used as a first guess to the velocity
        %field.
        disp('  -> Obtaining initial displacement field guess from PIV')
        if i<size(ims,2)
            pivout=usepivestimate(ims(i).cdata,ims(i+1).cdata,setfile);
            particles(i).pivguess=pivout;
        end
        %now we need to interpolate these velocities onto the particle
        %positions.
        disp(['  -> Interpolating velocities from PIV measurements onto', ...
            ' particle positions'])
        for jj=1:length(particles(i).blobs.area)
            particles(i).blobs.pivvel(jj,1)=interp2(particles(i).pivguess.x,...
                particles(i).pivguess.y,...
                particles(i).pivguess.fu,...
                particles(i).blobs.centr(jj,1),...
                particles(i).blobs.centr(jj,2));
            particles(i).blobs.pivvel(jj,2)=interp2(particles(i).pivguess.x,...
                particles(i).pivguess.y,...
                particles(i).pivguess.fv,...
                particles(i).blobs.centr(jj,1),...
                particles(i).blobs.centr(jj,2));
        end
    end
    
    %initialize the numofmatched variable
    numofmatched=0;
    
    if i>=2
        %disp(['  -> Interpolating velocities from  matched particles at' ...
        %	  ' previous timestep onto all particle positions'])
        
        fprintf('  -> Matching particles...')
        particles(i).alpha=matchparticles(particles(i-1), particles(i),setfile);
        numofmatched=sum(~isnan(particles(i).alpha));
        fprintf([num2str(numofmatched),' particles matched \n'])
        
        %Now we have to interpolate the velocities onto a regular grid,
        %as well as onto the new particle positions.
        alp=particles(i).alpha;
        nanarray=~isnan(alp);
        aind=find(~isnan(alp));
        disp('  - Interpolating gridded velocity field')
        
        uu=griddata(particles(i-1).blobs.centr(nanarray,1),...
            particles(i-1).blobs.centr(nanarray,2),...
            (particles(i).blobs.centr(alp(nanarray),1)-...
            particles(i-1).blobs.centr(nanarray,1))/dt,...
            particles(1).pivguess.x,particles(1).pivguess.y);
        vv=griddata(particles(i-1).blobs.centr(nanarray,1),...
            particles(i-1).blobs.centr(nanarray,2),...
            (particles(i).blobs.centr(alp(nanarray),2)-...
            particles(i-1).blobs.centr(nanarray,2))/dt,...
            particles(1).pivguess.x,particles(1).pivguess.y);
        [uu,vv]=naninterp(uu,vv);
        [uu,vv]=localfilt(particles(1).pivguess.x,particles(1).pivguess.y, ...
            uu,vv,3,'median');
        [uu,vv]=naninterp(uu,vv);
        particles(i-1).ptvgridvel.u=uu;
        particles(i-1).ptvgridvel.v=vv;
        %next we calculate the velocities of the matched particles
        particles(i-1).blobs.velocity(:,1)=(particles(i).blobs.centr(alp(aind),1)-...
            particles(i-1).blobs.centr(aind,1))/dt;
        particles(i-1).blobs.velocity(:,2)=(particles(i).blobs.centr(alp(aind),2)-...
            particles(i-1).blobs.centr(aind,2))/dt;
        
        %next we interpolate the gridded velocities of the tracked particles onto
        %the particle positions at frame i, to be used as an estimate for
        %the displacement at the next frame
        disp(['  - Generating displacement estimate for all particles for use' ...
            ' at next time step'])
        particles(i).blobs.ptvvel(:,1)=interp2(particles(1).pivguess.x,...
            particles(1).pivguess.y,...
            particles(i-1).ptvgridvel.u,...
            particles(i).blobs.centr(:,1),...
            particles(i).blobs.centr(:,2));
        particles(i).blobs.ptvvel(:,2)=interp2(particles(1).pivguess.x,...
            particles(1).pivguess.y,...
            particles(i-1).ptvgridvel.v,...
            particles(i).blobs.centr(:,1),...
            particles(i).blobs.centr(:,2));
    end
    
    % Here we display particle positions as we go
    if strcmp(showfigs,'yes')
        if i==1, %initialise the figure at the first frame
            figure, set(gcf,'Position',[80 20 600 700]);
            h1=subplot(2,2,1); colormap(hsv); title('Current image')
            %text will be output in this frame:
            h3=subplot(2,2,2); axis on; box on; title('Status')
            set(gca,'Color',[1 1 0.7]),set(gca,'XTick',[]),set(gca,'yTick',[])
            h2=subplot(2,1,2); set(gca,'Color',[0 0 0]), axis ij, hold on
        end
        subplot(2,2,1), imagesc(ims(i).cdata), hold on
        for ii=1:length(blobs.area)
            plot([blobs.bound(ii,1) blobs.bound(ii,1) blobs.bound(ii,1)+blobs.bound(ii,3)...
                blobs.bound(ii,1)+ blobs.bound(ii,3) blobs.bound(ii,1)],...
                [blobs.bound(ii,2) blobs.bound(ii,2)+blobs.bound(ii,4) blobs.bound(ii,2)+...
                blobs.bound(ii,4) blobs.bound(ii,2) blobs.bound(ii,2)],'k-');
        end
        %plot(blobs.centr(:,1),blobs.centr(:,2),'k.'), hold off,
        if i==1, myax=axis, else axis(myax), end %make sure the axis don't change
        
        subplot(2,2,2),
        if i>1, cla %delete(ht1,ht2,ht3,ht4);
        end
        ht1=text(0.02,0.95,['Frame number: ',num2str(i)]);
        ht2=text(0.02,0.85,['Number of particles: ',num2str(length(blobs.area))]);
        ht3=text(0.02,0.80,[' -Matched particles: ',num2str(numofmatched)]);
        ht4=text(0.02,0.75,[' -Matched 2 steps back: ']);
        subplot(2,1,2),
        h(comtel)=plot(blobs.centr(:,1),blobs.centr(:,2),'w.');
        drawnow, axis(myax)
        title('Visualisation of particle paths')
        if length(h)>1
            for j=length(h):-1:max([length(h)-comlen-1 1])
                set(h(j),'Color',[0 0 0]+(j-1)./min([comlen length(h)]))
            end
        end
        if length(h)==comlen
            delete(h(1))
            h(1:comlen-1)=h(2:comlen);
        else
            comtel=comtel+1;
        end
    end
end