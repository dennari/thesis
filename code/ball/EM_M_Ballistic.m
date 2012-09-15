function [new] = EM_M_Ballistic(p,m0T,gi,N,I1,I2,I3)
% em_lb_harmonic - lower bound and its gradient 
%
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement variance
% gi        a vector of indices to p, specifying the ones that are varying

global dt
    
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
          I = I2(1:3,1:3);
          new(3) = (6*I(1,1) - 3*I(1,2)*dt - 3*I(2,1)*dt + 2*I(2,2)*dt^2)/...
                   (N*dt^3);
        end
    end
    if(gi(j)==4) % dlb/dqy
          I = I2(3:4,3:4);
          new(4) = (6*I(1,1) - 3*I(1,2)*dt - 3*I(2,1)*dt + 2*I(2,2)*dt^2)/...
                   (N*dt^3);
    end
    if(gi(j)==5) % dlb/dr
      new(5) = (I3(1,1)+I3(2,2))/(2*N);
    end
end


end




