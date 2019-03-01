function [updatedblobs]=compareblobs(blobs,newblobs,dista,edgedista)
% COMPAREBLOBS - Check if blobs are unique in a grayscale image
% function [updatedblobs]=compareblobs(blobs,newblobs,maxdist)
%
% Checks if the newblobs are the same as the ones found at a lower
% threshold by checking if they are within a given distance from 
% each other


if nargin==2
    dista=3; % use distance=1 in test cases
    edgedista=1;
elseif nargin==3
  edgedista=1;
end

% The DISTA parameter is currently used as a measure of the distance
% between the centroid of the particles. Consider adding a measure of
% the overlap of the bounding box instead. This should be more
% accurate. Maybe a combination?


% Now we need to calculate the distance from all blobs to all newblobs.

oldx=blobs.centr(:,1);
oldy=blobs.centr(:,2);
newx=newblobs.centr(:,1)';
newy=newblobs.centr(:,2)';

%particle radius:
newr=sqrt(newblobs.area/pi);
oldr=sqrt(blobs.area/pi);
%

newindex=[];
for ii=1:length(newx)
  %Check the centroids
  dx=newx(ii)-oldx;   
  dy=newy(ii)-oldy; 
  dr=newr(ii)-oldr;
  dist=sqrt(dx.^2 + dy.^2); %distance between centroids
  
  if ~any(dist==0) & ~any(dist <=dista)
    angl1=asin((dy./dist)); %angle between centroids seen from oldblob
    angl2=acos(dx./dist);

    % now calculate the distance between "bounding boxes" (assumed
    % circular - but this is not entirely correct :-)
    %space=sqrt( (dx - (newr(ii).*cos(angl1)+oldr.*cos(angl1))).^2 +...
    %		(dy - (newr(ii).*sin(angl1)+oldr.*sin(angl1))).^2 );
    %
    % Therefor the spacex,spacey thing below which honors the
    % ellipticity of the blobs
    spacex=newx(ii)+newblobs.bound(ii,3)/2 - (oldx+blobs.bound(:,3)/2);
    spacey=newy(ii)+newblobs.bound(ii,4)/2 - (oldy+blobs.bound(:,4)/2);
    space=sqrt(spacex.^2 + spacey.^2);
    %[angl1, angl2, space, dist]
    
    % We also need to check the overlap of the bounding boxes. Perhaps
    % something ala "ismember" might work
    % any(any(ismember(newblobs.bound(ii),blobs.bound)))

    if ~any(space <= edgedista)  &...
	  ~any(newblobs.bound(ii,1) <= blobs.centr(:,1) &...
	  blobs.centr(:,1) <= newblobs.bound(ii,1)+newblobs.bound(ii,3) & ...
	  newblobs.bound(ii,2) <= blobs.centr(:,2) &...
	  blobs.centr(:,2) <= newblobs.bound(ii,2)+newblobs.bound(ii,4))
      %The last check above eliminates particles that have the
      %centroid of an older particle within their bounding box.
      
      newindex=[newindex,ii];
      %disp('Stays')
      %[angl1, angl2, space, dist]
    end 
  end
end

% The following matrix method uses a LOT of memory - therefore the loop above
%oldx2=repmat(oldx,1,length(newx));
%newx2=repmat(newx,length(oldx),1);
%oldy2=repmat(oldy,1,length(newy));
%newy2=repmat(newy,length(oldy),1);

%tmpx=oldx2 - newx2;
%tmpy=oldy2 - newy2;
%tmpm=sqrt(tmpx.^2 + tmpy.^2);

%newindex= ~any(tmpm <= dista);
%toc
updatedblobs.area=[blobs.area; newblobs.area(newindex)];
updatedblobs.centr=[blobs.centr; newblobs.centr(newindex,:)];
updatedblobs.bound=[blobs.bound; newblobs.bound(newindex,:)];

