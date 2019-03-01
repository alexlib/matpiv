function [x,y,u,v,scale]=ptvvekplot(p,field,time,scale,col);
% PTVVEKPLOT - plot vectors as arrows
%
% This function takes the output from MATPTV as an input-argument
% along with specification of which vector field to plot
%
% example -  this will plot the gridded velocity:
%
% >> ptvvekplot(particles,'ptvgridvel',2,0.01);
%
% and this will extract the velocity-field for t=2:
%
% >> [x,y,u,v]=ptvvekplot(particles,'ptvgridvel',2);
%

if nargin<3, disp('Wrong input to PTVVEKPLOT'); return, end

if nargin>=3, 
    if strcmp(field,'ptvgridvel')
        tmp=getfield(p(time),field);
        u=tmp.u; v=tmp.v;
        x=getfield(p(1).pivguess,'x'); y=getfield(p(1).pivguess,'y');
    elseif strcmp(field,'pivguess')
        tmp=getfield(p(time),field);
        u=tmp.fu; v=tmp.fv;
        x=getfield(p(1).pivguess,'x'); y=getfield(p(1).pivguess,'y');
    elseif strcmp(field,'velocity') | strcmp(field,'ptvvel')
        tmp=getfield(p(time).blobs,field);
        u=tmp(:,1); v=tmp(:,2);
        
        tmp=getfield(p(time).blobs,'centr');
        alp=getfield(p(time+1),'alpha');
        if strcmp(field,'velocity')
            x=tmp(~isnan(alp),1);
            y=tmp(~isnan(alp),2);
        else
            x=tmp(:,1); y=tmp(:,2);
        end
    else
       disp('Wrong input to PTVVEKPLOT -  no such velocity field'); return 
    end
  if nargin==3 
    scale=3/max(sqrt(u(:).^2 + v(:).^2)); col='r';
  elseif (nargin==4 & isstr(scale))
    col=scale; scale=3/max(sqrt(u(:).^2 + v(:).^2));
  end
  
end

if nargout==0
  vekplot2(x/10,y/10,u,v,scale,col);
elseif nargout==1
  vekplot2(x/10,y/10,u,v,scale,col);
  x=scale;
elseif nargout>0 & nargout<4
  disp('Wrong number of output-arguments')
end