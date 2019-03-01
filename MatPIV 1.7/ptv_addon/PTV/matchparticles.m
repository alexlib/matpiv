function [alpha]=matchparticles(b1,b2,varargin)
% MATCHPARTICLES - match particles in MATPTV
%
%
%
%
if nargin==2
  m=2; %might introduce this in the input arguments
  maxerror=2; %max squared error in position estimate
else
  setfile=deal(varargin{:});
  eval(setfile)
end


if isfield(b1.blobs,'ptvvel')
  if ~isempty(b1.blobs.ptvvel)
  %  if isempty(b1.blobs.pivvel)
      ui=b1.blobs.ptvvel(:,1); vi=b1.blobs.ptvvel(:,2);
  %  else
  %    ui=whodoyoutrust(1)*b1.blobs.ptvvel(:,1) + ...
%	 whodoyoutrust(2)*b1.blobs.pivvel(:,1); 
%      vi=whodoyoutrust(1)*b1.blobs.ptvvel(:,2) + ...
%	 whodoyoutrust(2)*b1.blobs.pivvel(:,2);
%    end
  else
    ui=b1.blobs.pivvel(:,1); vi=b1.blobs.pivvel(:,2);
  end
else
  ui=b1.blobs.pivvel(:,1); vi=b1.blobs.pivvel(:,2);
end

xi=b1.blobs.centr(:,1); yi=b1.blobs.centr(:,2);
xj=b2.blobs.centr(:,1); yj=b2.blobs.centr(:,2);

dt=b2.t-b1.t;
Bij=zeros(length(xj),length(xi));
% first approach: match with whichever particle is closest in frame 2
for i=1:length(xi)
  for j=1:length(xj)
    Bij(j,i) = (sqrt( (xi(i) + ui(i)*dt - xj(j)).^2 + (yi(i) + vi(i)*dt - yj(j)).^2 ) ).^m; 
  end
  if ~any(isnan(Bij(:,i)))
    bminind=find(Bij(:,i)==min(Bij(:,i)));
    %if length(bminind)>1
      % this means we have more than one matching particle in frame 2
      % here we'll match using particle properties, not just displacement
      
    %end
    if Bij(bminind,i)<maxerror
      alpha(i)=bminind;
    else
      alpha(i)=NaN;
    end
  end
end

%Now check if some particles are matched with more than one in frame
%2
alpha(alpha==0)=nan; %no indexes should be 0
for i=1:length(alpha)
  [tind]=find(alpha==alpha(i));
  if length(tind)>1
    b_ind=find(Bij(alpha(tind),tind)==min(min(Bij(alpha(tind),tind)))); 
    if length(b_ind)>1
      alpha(b_ind)=NaN;
    else 
      alpha(tind(~b_ind))=NaN;
    end
  end
end
%minimization
%her er det noe muffens med dimensjonene.
%alpha=b1.blobs.centr'*b2.blobs.centr \ Bij;

