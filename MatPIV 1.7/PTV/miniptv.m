
minstepart=4; % minste partikkel vi tillater
storpart=40; %største partikkel vi tillater før den forkastes
maxavstand=5; % max avstand (i piksler) før vi forkaster matchen
              % hvis hastighetene er store og vi får få treff, må
              % denne paremeteren økes. 
p1=blobrec(a,graythresh(a),minstepart,storpart);
p2=blobrec(b,graythresh(a),minstepart,storpart);
xi=p1.centr(:,1); yi=p1.centr(:,2);
xj=p2.centr(:,1); yj=p2.centr(:,2);
Bij=zeros(length(xj),length(xi));
for i=1:length(xi)
  for j=1:length(xj)
Bij(j,i)=(sqrt((xi(i)-xj(j)).^2 + (yi(i)-yj(j)).^2 )).^2;
 end
  if ~any(isnan(Bij(:,i)))
    bminind=find(Bij(:,i)==min(Bij(:,i)));
    if Bij(bminind,i)<maxavstand^2
      alpha(i)=bminind;
    else
      alpha(i)=NaN;
    end
  end
end
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

imagesc(abs(a-b))
colormap(gray)
hold on
isa=isnan(alpha);
x=p1.centr(~isa,1); y=p1.centr(~isa,2);
u=p2.centr(alpha(~isa),1)-p1.centr(~isa,1);
v=p2.centr(alpha(~isa),2)-p1.centr(~isa,2);
vekplot2(x,y,u,v,1)