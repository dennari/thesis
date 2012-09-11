function [new] = EM_M_Ballistic(p,m0T,gi,N,I1,I2,I3)
% em_lb_harmonic - lower bound and its gradient 
%
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement variance
% gi        a vector of indices to p, specifying the ones that are varying

%global P0 g0x g0y
    
%Q = ballisticQ(p(3),p(4));
%R = p(5)^2*eye(size(I3,1));
%m0 = [0 p(1) g0x 0 p(2) g0y]';

% gradient
new = p;
for j=1:numel(gi)
    % assume all zero
    
    if(gi(j)==1) % dlb/dqw
    end
    if(gi(j)==3) % dlb/dqx
        if iscell(p)
          new{1} = Q;
        else
          dQ = ballisticQ(1,0); % the multipliers
          qI = 1/N*I2(1:3,1:3)/dQ(1:3,1:3);
          diag(qI)
          new(3) = qI(3,3);
        end
    end
    if(gi(j)==4) % dlb/dqy
          dQ = ballisticQ(0,1); % the multipliers
          t = 1/N*I2(4:6,4:6)/dQ(4:6,4:6);
          new = abs(t(3,3));
    end
    if(gi(j)==5) % dlb/dr
    end
end


end




