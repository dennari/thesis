function [lb,glb] = em_lb_harmonic(p,vi,Y,JM,JS)
% em_lb_harmonic - lower bound and its gradient 
%
% p{1}=qw   angular velocity variance
% p{2}=qx   OR Qx \in 2x2, the signal component variance 
% p{3}=r    measurement variance
% p{4}=m0   prior mean
% p{5}=P0   prior covariance matrix
% p{6}=H    state-to-measurement matrix, fixed
% p{7}=dt   the discretization time interval
% vi        a vector of indices to p, specifying the ones that are varying
% Y         the measurements
% JM        cross timestep smoothed means
% JP        cross timestep smoothed covariance matrices

[qw,qx,R,m0,P0,H,dt] = p{:};

N = size(Y,2);
xDim = size(m0,1);
yDim = size(Y,1);
c = (xDim-1)/2; % number of harmonics
if numel(qx) == 1 && c > 1
    qqx = qx*ones(1,c);
else
    qqx = qx;
end
Q = sinusoid_Q(qw,qqx,dt);



[eij,wij] = CKFPoints(2*xDim);
M = size(eij,2);
% weight matrix
wmx = repmat(wij(1,:),xDim,1);
wmy = repmat(wij(1,:),yDim,1);

I1 = zeros(xDim,xDim);
I2 = zeros(xDim,xDim);
I3 = 0;	

%newer = 1:xDim;
%older = (xDim+1):(2*xDim);
older = 1:xDim;
newer = (xDim+1):(2*xDim);


%IIx = [eye(xDim);-eye(xDim)];
%IIy = [eye(yDim);-eye(yDim)];

for k=1:N
    % compute integrals I_2 and I_3 that are w.r.t N(m3,P3)
	m = JM(:,k);
	S = JS(:,:,k);
    
    if k == 1
        % p(x_0|y_1:N)
        d = m(older)-m0;
        I1 = S(older,older)+d*d';
    end
    
    % get the sigma points from the joint distribution
    [xk1,xk] = jsig(m,S);
    
	d = xk-sinusoid_f(xk1,dt);
	I2 = I2 + (wmx.*d)*d';
	
    %j = [xk ; model.f(xk1,k,p)];
	%I2 = I2 + IIx'*(wmx.*j)*j'*IIx;
	
	d = repmat(Y(:,k),1,M)-H*xk;
    I3 = I3 + (wmy.*d)*d';

end


% gradient

% in the case that numel(qx) == num of components,
% we will put each dlb/dqx(i)
% one after the other in glb
% otherwise we will have just dlb/dqx(1)


% expand the "2" in vi to c*"2"
if sum(vi==2) == 1 && numel(qx) == c
    glb = zeros(numel(vi)+c-1,1);
    vii = glb;
    qi = find(vi==2);
    if (qi > 1)
        vii(1:(qi-1)) = vi(1:(qi-1));
    end
    vii(qi:(qi+c-1)) = 2;
    if qi < numel(vi) 
        vii((qi+c):end) = vi(qi:end);
    end
    dqxi = eye(c);
else % the qx variances are the same for all components
    glb = zeros(numel(vi),1);
    vii = vi;
    dqxi = ones(1,c);
end

ri = 1; % another counter
for j=1:numel(glb)
    % assume all zero
    dSig = zeros(size(P0)); dQ = zeros(size(Q));dR = zeros(size(R)); dmu = zeros(xDim,1);
    
    if(vii(j)==1) % dlb/dqw
        dQ = zeros(size(Q));
        dQ(1,1) = 1;
    end
    if(vii(j)==2) % dlb/dqx(ri)
        dQ = sinusoid_Q(0,dqxi(ri,:),dt);
        ri = ri + 1;
    end
    if(vii(j)==3) % dlb/dr
        dR = 1;
    end
    
    % x_0
    glb1 = P0\(dSig/P0*I1+2*dmu*(JM(older,1)-m0)'+dSig);
    % x_k|x_(k-1)
    glb2 = Q\(dQ/Q*I2+N*dQ);
    % y_k|x_k
    glb3 = R\(dR/R*I3+N*dR);

    glb(j) = 0.5*(trace(glb1)+trace(glb2)+trace(glb3));
end


I1 =  log(det(P0))+trace(I1/P0);
I2 = N*log(det(Q))+trace(I2/Q);
I3 = N*log(R)+trace(I3/R);

lb = 0.5*(I1+I2+I3);


function [old,new]=jsig(m,S)
    mr = repmat(m,1,M);
    xi = mr+chol(S,'lower')*eij;
	% sigma points for x_k
    new = xi(newer,:);
    % sigma points for x_(k-1)
	old = xi(older,:);
end

end




