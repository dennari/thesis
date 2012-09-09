function [lb,glb] = EM_LB_Ballistic(p,vi,Y,JM,JP)
% em_lb_harmonic - lower bound and its gradient 
%
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement variance
% vi        a vector of indices to p, specifying the ones that are varying
% Y         the measurements
% JM        cross timestep smoothed means
% JP        cross timestep smoothed covariance matrices
global P0 H A g0x g0y

[v0x,v0y,qx,qy,R] = p{:};

N = size(Y,2);
xDim = size(P0,1);
yDim = size(Y,1);
m0 = [0 v0x g0x 0 v0y g0y]';
Q = ballisticQ(qx,qy);

I1 = 0;
I2 = 0;
I3 = 0;	

%newer = 1:xDim;
%older = (xDim+1):(2*xDim);
older = 1:xDim;
newer = (xDim+1):(2*xDim);


for k=1:N
    % compute integrals I_2 and I_3 that are w.r.t N(m3,P3)
	m = JM(:,k);
	P = JP(:,:,k);
    
    
    if k == 1
        % p(x_0|y_1:N)
        d = m(older)-m0;
        I1 = P(older,older)+d*d';
    end
    
    I2 = I2 + P+m*m';
    
    y = Y(:,k);
    I3 = I3 + [y*y'         y*m(newer)';
               m(newer)*y'  m(newer)*m(newer)'];

end

IA = [eye(xDim) -A]';
IH = [eye(yDim) -H]';
I2 = IA'*I2*IA;
I3 = IH'*I3*IH;

% gradient

glb = zeros(numel(vi));
for j=1:numel(vi)
    % assume all zero
    dSig = zeros(size(P0)); dQ = zeros(size(Q));dR = zeros(size(R));
    dmu = zeros(xDim,1);
    
    if(vi(j)==1) % dlb/dqw
        dQ = zeros(size(Q));
        dQ(1,1) = 1;
    end
    if(vi(j)==3) % dlb/dqx
        dQ = ballisticQ(1,0);
    end
    if(vi(j)==4) % dlb/dqy
        dQ = ballisticQ(0,1);
    end
    if(vi(j)==5) % dlb/dr
        dR = eye(yDim);
    end
    
    % x_0
    glb1 = P0\(dSig/P0*I1+2*dmu*(JM(older,1)-m0)'+dSig);
    % x_k|x_(k-1)
    glb2 = Q\(dQ/Q*I2+N*dQ);
    % y_k|x_k
    glb3 = R\(dR/R*I3+N*dR);

    glb(j) = 0.5*(trace(glb1)+trace(glb2)+trace(glb3));
end




I1 =   log(det(P0))+trace(I1/P0);
I2 = N*log(det(Q))+trace(I2/Q);
I3 = N*log(det(R))+trace(I3/R);

lb = -0.5*(I1+I2+I3);

end




