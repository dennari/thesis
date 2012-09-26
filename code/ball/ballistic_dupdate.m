function [dm,dP,dmy,dSy ] = ballistic_dupdate(gi,dm_,dP_,P_,K,my,Sy,y,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
			% for AR(1) model

    global H

    C = P_*H';

    dR = zeros(size(H,1));

    if(gi==3) % dlh/dr
        dR = eye(size(dR,1))*2*exp(2*p(3));
    end
    
    
    
    dmy = H*dm_;
    dSy = H*dP_*H' + dR;
    dK = dP_*H'/Sy - C*(Sy\dSy/Sy);
    
    
    dm = dm_+dK*(y-my)-K*dmy; 
    dP = dP_-dK*Sy*K'-K*dSy*K'-K*Sy*dK'; 

    %do0 = trace(Sy\dSy);
    %do1 = (H*dm_)'/Sy*d;
    %do2 = d'*(Sy\dSy/Sy)*d;

    % gradient
    %glh(j) =  (0.5*(do2-do0)+do1);

end
    

            


