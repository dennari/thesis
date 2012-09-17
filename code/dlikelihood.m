function [dlh] = dlikelihood(Sy,dSy,y,my,dmy)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
			% for AR(1) model
 
    d = y-my;
    
    do0 = trace(Sy\dSy);
    do1 = dmy'/Sy*d;
    do2 = d'*(Sy\dSy/Sy)*d;

    % partial derivative
    dlh = (0.5*(do2-do0)+do1);

end
    

            


