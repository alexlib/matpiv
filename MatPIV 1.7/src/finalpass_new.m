function [xp,yp,up,vp,SnR,Pkh,brc]=finalpass_new(A,B,N,ol,idx,idy,Dt,maske)
% function [x,y,u,v,SnR,PeakHeight,brc]=finalpass_new(A,B,N,ol,idx,idy,Dt,mask)
%
% Provides the final pass to get the displacements with
% subpixel resolution. Uses sub-pixel displacements in Fourier Space
% following the paper by Qian and Cowen (2005, Experiments in Fluids).
%
%

% 1999 - 2011, J. Kristian Sveen (jks@math.uio.no)
% For use with MatPIV 1.7a, Copyright
% Distributed under the terms of the GNU - GPL license
% timestamp: 09:40, 4 Mar 2011

if length(N)==1
    M=N;
else
    M=N(1); N=N(2);
end
cj=1;
[sy,sx]=size(A);

% Allocate space for matrixes
xp=zeros(ceil((size(A,1)-N)/((1-ol)*N))+1, ...
    ceil((size(A,2)-M)/((1-ol)*M))+1);
yp=xp; up=xp; vp=xp; brc=xp; SnR=xp; Pkh=xp;

% Variables used for subpixel displacement in Fourier Domain
min_res=0.005; % minimum residual to reach before breaking out of
% sub-pixel iterations
breakoutcounter=1; % count the iterations
max_iterations=10; % max iterations before breaking out
I=1:2*M; I(I>M)=I(I>M)-2*M; I=repmat(I,2*N,1); % used in the sub-pixel
% window shift
J=(1:2*N)'; J(J>N)=J(J>N)-2*N; J=repmat(J,1,2*M);
W=weight('cosn',[M,N],20); % weights used in the sub-pixel window shift
%W2=1-weight('cosn',2*[M N],20); %weight for FFT filtering

if nargin==8, 
    if ~isempty(maske)
        IN=zeros(size(maske(1).msk));
        for ii=1:length(maske)
            IN=IN+double(maske(ii).msk);
        end
    else IN=zeros(size(A)); 
    end, 
end

fprintf([' Continuous windows shifting in Fourier space\n', ...
    '  Local iterations applied\n',...
    '  - Using ',num2str(M),'*',num2str(N),...
    ' interrogation windows! \n'])
%%%%%%%%%%%%%%% MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
tic
for jj=1:((1-ol)*N):sy-N+1
    ci=1;
    for ii=1:((1-ol)*M):sx-M+1
        if IN(jj+N/2,ii+M/2)~=1
            if isnan(idx(cj,ci))
                idx(cj,ci)=0;
            end
            if isnan(idy(cj,ci))
                idy(cj,ci)=0;
            end
            if jj+idy(cj,ci)<1
                idy(cj,ci)=1-jj;
            elseif jj+idy(cj,ci)>sy-N+1
                idy(cj,ci)=sy-N+1-jj;
            end
            if ii+idx(cj,ci)<1
                idx(cj,ci)=1-ii;
            elseif ii+idx(cj,ci)>sx-M+1
                idx(cj,ci)=sx-M+1-ii;
            end
            D2=B(jj+idy(cj,ci):jj+N-1+idy(cj,ci),ii+idx(cj,ci):ii+M-1+idx(cj,ci));
            E=A(jj:jj+N-1,ii:ii+M-1);
            stad1=std(E(:));
            stad2=std(D2(:));
            if stad1==0, stad1=1; end
            if stad2==0, stad2=1; end
            
            % use weights
            E=(E-mean(E(:))).*W;
            F=(D2-mean(D2(:))).*W;
            E=E-mean(E(:));
            F=F-mean(F(:));
            
            % take zero-padded Fourier Transform
            mf = 2^nextpow2(M+N);
            nf = mf;
            at = fft2(E,nf,mf);
            bt = fft2(conj(F(end:-1:1,end:-1:1)),nf,mf);
            %no zero-padding version - need to change I and J if this
            %is uncommented
            %at = fft2(E);
            %bt = fft2(F);
            %%%%%%%%%%%%%%%%%%%%%% Calculate the normalized correlation:
            R=real(ifft2(bt.*at));
            R(end,:)=[]; R(:,end)=[];
            R=real(R)./(N*M*stad1*stad2);
            %%%%%%%%%%%%%%%%%%%%%% Find the position of the maximal value of R
            %%%%%%%%%%%%%%%%%%%%%% _IF_ the standard deviation is NOT NaN.
            if all(~isnan(R(:))) && ~all(R(:)==0)  %~isnan(stad1) & ~isnan(stad2)
                if size(R,1)==(N-1)
                    [max_y1,max_x1]=find(R==max(R(:)));
                    
                else
                    [max_y1,max_x1]=find(R==max(max(R(0.5*N+2:1.5*N-3,...
                        0.5*M+2:1.5*M-3))));
                end
                if length(max_x1)>1
                    max_x1=round(sum(max_x1.^2)./sum(max_x1));
                    max_y1=round(sum(max_y1.^2)./sum(max_y1));
                end
                
                % loop on integer basis to make sure we've converged to
                % +-0.5 pixels before entering subpixel shifting
                stx=idx(cj,ci); sty=idy(cj,ci);
                
                while max_x1~=M && max_y1~=N && ...
                        breakoutcounter<max_iterations &&...
                        jj+sty>0 && ii+stx>0 && ii+M-1+stx<=sx && jj+N-1+sty<=sy
                    D2=B(jj+sty:jj+N-1+sty,...
                        ii+stx:ii+M-1+stx);
                    F=(D2-mean(D2(:))).*W; F=F-mean(F(:));
                    bt = fft2(conj(F(end:-1:1,end:-1:1)),nf,mf);
                    R=ifft2(bt.*at);
                    R(end,:)=[]; R(:,end)=[];
                    R=real(R)./(N*M*stad1*stad2);
                    [max_y1,max_x1]=find(R==max(R(:)));
                    stx=stx + (M-max_x1);
                    sty=sty + (N-max_y1);
                    
                    breakoutcounter=breakoutcounter+1;
                end
                
                if breakoutcounter~=max_iterations
                    %update these only IF convergence was met, that is, we
                    %used less than max_iterations
                    idx(cj,ci)=stx; idy(cj,ci)=sty; 
                    breakoutcounter=1; % only reset if converged
                end
                
                %Only enter next bit if the peak is not located at the
                %edges of the correlation plane
                if max_x1~=1 && max_y1~=1 && max_x1~=M-1 && max_y1~=N-1
                    % 3-point peak fit using centroid, gaussian (default)
                    % or parabolic fit
                    [x0 y0]=intpeak(max_x1,max_y1,R(max_y1,max_x1),...
                        R(max_y1,max_x1-1),R(max_y1,max_x1+1),...
                        R(max_y1-1,max_x1),R(max_y1+1,max_x1),2,[M,N]);
                    X0=x0; Y0=y0;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % here we do the subpixel shifts in Fourier space
                    while (abs(x0)>min_res || abs(y0)>min_res) &&...
                            breakoutcounter<max_iterations
                        bt2=(exp(2*1i*pi*( (I*(X0)/size(I,2)) + ...
                            (J*(Y0)/size(J,1))))).*bt;
                        % At this point we could do some simple
                        % FFT-filtering like Todd and Liao suggests in
                        % their paper
                        % for example only keeping a certain number of
                        % Fourier components
                        % 26/3-05: This makes the peak a lot wider/broader:
                        % W2=1-weight('cosn',64,20)
                        % R=ifft2( bt2.*at.*W2 );
                        %
                        R=ifft2( bt2.*at );
                        R(end,:)=[]; R(:,end)=[];
                        R=real(R)./(N*M*stad1*stad2);
                        
                        [dy,dx]=find(R==max(R(:)));
                        X0=X0+(M-dx); Y0=Y0+(N -dy);
                        if dx>1 && dx<2*M-1 && dy>1 && dy<2*N-1
                            %only gaussian fit here
                            x0= -(log(R(dy,dx-1))-log(R(dy,dx+1)))/...
                                (2*log(R(dy,dx-1))-4*log(R(dy,dx))+2*log(R(dy,dx+1)));
                            y0= -(log(R(dy-1,dx))-log(R(dy+1,dx)))/...
                                (2*log(R(dy-1,dx))-4*log(R(dy,dx))+2*log(R(dy+1,dx)));
                            
                            X0=X0-x0; Y0=Y0-y0;
                            breakoutcounter=breakoutcounter+1;
                            %if ~isreal(x0) | ~isreal(y0)
                            %  disp([num2str([ci,cj]),', ',...
                            %	    num2str([stad1 stad2 std(A(:))])])
                            %    end
                        else
                            X0=nan; Y0=nan;
                            breakoutcounter=16;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
                    % Find the signal to Noise ratio
                    R2=R;
                    try
                        R2(max_y1-3:max_y1+3,max_x1-3:max_x1+3)=NaN;
                    catch
                        R2(max_y1-1:max_y1+1,max_x1-1:max_x1+1)=NaN;
                    end
                    if size(R,1)==(N-1)
                        [p2_y2,p2_x2]=find(R2==max(R2(:)));                        
                    else
                        [p2_y2,p2_x2]=find(R2==max(max(R2(0.5*N:1.5*N-1,0.5*M:1.5*M-1))));
                    end
                    if length(p2_x2)>1
                        p2_x2=p2_x2(round(length(p2_x2)/2));
                        p2_y2=p2_y2(round(length(p2_y2)/2));
                    elseif isempty(p2_x2)
                        
                    end
                    % signal to noise:
                    snr=R(max_y1,max_x1)/R2(p2_y2,p2_x2);
                    % signal to mean:
                    %snr=R(max_y1,max_x1)/mean(R(:));
                    % signal to median:
                    %snr=R(max_y1,max_x1)/median(median(R(0.5*N+2:1.5*N-3,...
                    %    0.5*M+2:1.5*M-3)));
                    
                    %%%%%%%%%%%%%%%%%%%%%% Store the displacements, SnR and Peak Height.
                    up(cj,ci)=(-X0+idx(cj,ci))/Dt;
                    vp(cj,ci)=(-Y0+idy(cj,ci))/Dt;
                    xp(cj,ci)=(ii+(M/2)-1);
                    yp(cj,ci)=(jj+(N/2)-1);
                    SnR(cj,ci)=snr;
                    Pkh(cj,ci)=R(max_y1,max_x1);
                else
                    up(cj,ci)=NaN; vp(cj,ci)=NaN; SnR(cj,ci)=NaN; Pkh(cj,ci)=0;
                    xp(cj,ci)=(ii+(M/2)-1);
                    yp(cj,ci)=(jj+(N/2)-1);
                end
            else
                up(cj,ci)=NaN; vp(cj,ci)=NaN; SnR(cj,ci)=NaN; Pkh(cj,ci)=0;
                xp(cj,ci)=(ii+(M/2)-1);
                yp(cj,ci)=(jj+(N/2)-1);
            end
            ci=ci+1;
        else
            xp(cj,ci)=(M/2)+ii-1;
            yp(cj,ci)=(N/2)+jj-1;
            up(cj,ci)=NaN; vp(cj,ci)=NaN;
            SnR(cj,ci)=NaN; Pkh(cj,ci)=NaN;ci=ci+1;
        end
        brc(cj,ci)=breakoutcounter;
    end
    breakoutcounter=1;
    
    % disp([num2str((cj-1)*(ci)+ci-1) ' vectors in ' num2str(toc) ' seconds'])
    fprintf('\r No. of vectors: %d', ((cj-1)*(ci)+ci-1) -sum(isnan(up(:))))
    fprintf(', Seconds taken: %f', toc);
    cj=cj+1;
end

fprintf('\n')