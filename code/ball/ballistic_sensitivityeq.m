function [dm,dP,glh ] = ballistic_sensitivityeq(dm,dP,A,H,gi,d,S,C,m,P,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
			% for AR(1) model

K = C/S;

glh = zeros(numel(gi),1);            
    
for j=1:numel(gi)
    dA = zeros(size(A));
    dQ = zeros(size(A));
    dR = zeros(size(d,1));

    if(gi(j)==1) % dlh/dqx
        dQ = ballisticQ2D(1,0,1)*2*exp(2*p(1));
    end
    if(gi(j)==2) % dlh/dqy
        dQ = ballisticQ2D(0,1,1)*2*exp(2*p(2));
    end
    if(gi(j)==3) % dlh/dr
        dR = eye(size(dR,1))*2*exp(2*p(3));
    end
    
    

    dm_ = dA*m+A*dm;
    dP_ = dA*P*A'+A*dP*A'+A*P*dA'+dQ;
    
    dS = H*dP_*H' + dR;
    dK = dP_*H'/S - C*(S\dS/S);

    dm = dm_+dK*d-K*H*dm_; 
    dP = dP_-dK*S*K'-K*dS*K'-K*S*dK'; 

    do0 = trace(S\dS);
    do1 = (H*dm_)'/S*d;
    do2 = d'*(S\dS/S)*d;

    % gradient
    glh(j) =  (0.5*(do2-do0)+do1);

end
    

            


