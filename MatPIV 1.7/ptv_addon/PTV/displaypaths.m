function []=displaypaths(p,lenofpaths)
%
%
%
%

if nargin==1
  lenofpaths=3;
end
  
len=length(p);

if rem(lenofpaths,2)==0
  start=lenofpaths/2; stopp=len-lenofpaths/2-1;
else
  start=ceil(lenofpaths/2); stopp=len-ceil(lenofpaths/2)-1;
end
%[start stopp]

for i=start:stopp
  
  vecx=[]; vecy=[];aind=[];alp=[]; aind2=[];
  for j=1:lenofpaths
    alp1=p(i+j-1).alpha;
    alpind1=~isnan(alp1);
    aind1=find(alpind1);
    if j>1
      aind=intersect(aind1(:),aind);
    else 
      aind=aind1(:);
    end
    %alp=[alp;alp1];
  end   
  
    
  for j=1:lenofpaths   
    if j==1, 
      pp=p(i+j-1).blobs.centr(aind,1)'; 
      rr=p(i+j-1).blobs.centr(aind,2)'; 
    else, 
      pp=p(i+j-1).blobs.centr(p(i+j-1).alpha(aind),1)'; 
      rr=p(i+j-1).blobs.centr(p(i+j-1).alpha(aind),2)'; 
    end
    vecx=[vecx;pp, nan];
    vecy=[vecy;rr, nan];
  end
 
  plot(vecx,vecy,'k.-')
  
  %plot([p(i-1).blobs.centr(aind,1)',nan;p(i).blobs.centr(alp1(aind),1)',nan;...
  %	p(i+1).blobs.centr(alp2(aind),1)',nan;nan*ones(length(aind),1)',nan],...
  %       [p(i-1).blobs.centr(aind,2)',nan;p(i).blobs.centr(alp1(aind),2)',nan;...
  %	p(i+1).blobs.centr(alp2(aind),2)',nan;nan*ones(length(aind),1)',nan],'k.-')
  
  drawnow
    
end