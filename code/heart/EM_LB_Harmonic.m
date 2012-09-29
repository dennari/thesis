function [glb] = EM_LB_Harmonic(p,m0T,gi,N,I1,I2,I3)
% EM_LB_Harmonic - score function evaluation 
%
% parameters are 
% p(1)=lqw,    log angular velocity variance
% p(2)=lr,     log measurement variance
% p(3:3+c-1)   log component variances





% gradient

% in the case that numel(qx) == num of components,
% we will put each dlb/dqx(i)
% one after the other in glb
% otherwise we will have just dlb/dqx(1)


% expand the "3" in gi to c*"3"
% if sum(gi==3) == 1 && numel(p) == c+2
%     glb = zeros(numel(vi)+c-1,1);
%     vii = glb;
%     qi = find(vi==2);
%     if (qi > 1)
%         vii(1:(qi-1)) = vi(1:(qi-1));
%     end
%     vii(qi:(qi+c-1)) = 2;
%     if qi < numel(vi) 
%         vii((qi+c):end) = vi(qi:end);
%     end
%     dqxi = eye(c);
% else % the qx variances are the same for all components
%     glb = zeros(numel(vi),1);
%     vii = vi;
%     dqxi = ones(1,c);
% end

%ri = 1; % another counter
glb = zeros(numel(gi));
for j=1:numel(gi)

    
    if(gi(j)==1) % dlb/dqw
        %dQ(1,1) = 1*2*exp(2*p(gi(j)));
        qw = exp(2*p(1));
        glb(j) = I2(1,1)/qw-N;
        %Q = sinusoid_Q(p(1),p(3:end));
        %dQ = zeros(size(Q)); dQ(1,1) = 1;
        
        %glb2 = 0.5*trace((Q\(dQ/Q*I2-N*dQ)))*2*qw;
    end
    if(gi(j) >= 3) % dlb/dqx(ri)
        
        %ri = ri + 1;
        %wh = ones(1,c);
        %if sum(gi>=3) > 1
        %  wh = zeros(1,c);
        %  wh(gi-2) = 1;
        %end 
        qx = exp(2*p(3));
        dQ = sinusoid_Q(0,1,1);
        Q = sinusoid_Q(p(1),p(3));
        Q(2,2) = 1; Q(4,4) = 1; Q(6,6) = 1;
        glb(j) = trace(Q\(dQ/Q*I2-N*dQ))*qx;
    end
    if(gi(j)==2) % dlb/dr
        r = exp(2*p(2));
        glb(j) = I3/r-N;
    end
    
    % x_0
    %glb1 = P0\(dSig/P0*I1+2*dmu*(m0T-m0)'-dSig);
    % x_k|x_(k-1)
    %glb2 = Q\(dQ/Q*I2-N*dQ);
    % y_k|x_k
    %glb3 = R\(dR/R*I3-N*dR);

    %glb(j) = 0.5*(trace(glb2)+trace(glb3));
end


end




