function [dm_,dP_] = ballistic_dpred(gi,dm,dP,m,P,p)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
			% for AR(1) model

global dt A



    dA = zeros(size(A));
    dQ = zeros(size(A));
    du = zeros(size(A,1),1);

    if(gi==1) % dlh/dqx
        dQ = ballisticQ2D(1,0,1)*2*exp(2*p(1));
    end
    if(gi==2) % dlh/dqy
        dQ = ballisticQ2D(0,1,1)*2*exp(2*p(2));
    end
    if(gi==4) % dlh/dq (joint process noise)
        dQ = ballisticQ2D(1,1,1)*2*p(1);
    end
    if(gi==5 || gi==6) % dlb/du
        du(2*gi-8) = dt;        
    end
    

    dm_ = dA*m+A*dm+du;
    dP_ = dA*P*A'+A*dP*A'+A*P*dA'+dQ;
    


    

            


