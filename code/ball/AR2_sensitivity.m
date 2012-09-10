function [dm,dP,glh ] = AR2_sensitivity(dm,dP,A,H,vi,d,S,C,m,P)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
			% for AR(1) model

K = C/S;

glh = zeros(numel(vi),1);            
    
for j=1:numel(vi)

    dA = zeros(size(A));
    dQ = zeros(size(A));
    dR = zeros(size(d,1));

    if(vi(j)==1) % dlh/da1
        dA(1,1) = 1;
    end
    if(vi(j)==2) % dlh/da2
        dA(1,2) = 1;
    end
    if(vi(j)==3) % dlh/dq1
        dQ(1,1) = 1;
    end
    if(vi(j)==4) % dlh/dq2
        dQ(2,2) = 1;
    end
    if(vi(j)==5) % dlh/dr
        dR(1,1) = 1;
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
    

            


