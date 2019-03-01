load snr_ims

N=64;
tel2=1;
maxlen=length(1:48:size(A,2)-N);
tic
for ii=1:48:size(A,2)-N
  for jj=1:48:size(A,1)-N
    a=A(jj:jj+N-1,ii:ii+N-1);
    b=B(jj:jj+N-1,ii:ii+N-1);
    a=a-mean(a(:));
    b=b-mean(b(:));
    
    tel=1;  

    for i=0:0.01:1
      R=pcorr2(a,b,'nopad',i,'orig');
      [y1,x1]=find(R==max(R(:)));
      R2=R;
      R2(y1-3:y1+3,x1-3:x1+3)=NaN;
      [p2_y2,p2_x2]=find(R2==max(R2(:)));
      
      SnR(tel2,tel)=R(y1,x1)/R2(p2_y2,p2_x2);
      tel=tel+1;
    end
    toc
    stad(tel2)=std(a(:));%./std(A(:));
    skew(tel2)=skewness(a(:));
    kurt(tel2)=kurtosis(a(:));
    
    %if stad(tel2)>1, stad(tel2)=stad(tel2)-1; end
    %stad(tel2)=1-stad(tel2);
    
    tel2=tel2+1;
  end
  %fprintf([num2str(tel2/(N*maxlen)),', '])
end
fprintf('\n')

figure
plot([0:0.01:1],mean(SnR),'r+-')
hold on
plot([0:0.01:1],mean(SnR)-std(SnR),'.:')
plot([0:0.01:1],mean(SnR)+std(SnR),'.:')
plot([0.4250 0.4250],[0 9],'g-')
title('Signal to Noise ratio')
ylabel('SnR')
xlabel('\alpha')

legend('mean','mean-std','mean+std','mean \alpha used')
