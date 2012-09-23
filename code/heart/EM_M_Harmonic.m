function [new] = EM_M_Harmonic(p,m0T,gi,N,I1,I2,I3)
% EM_M_Harmonic - Maximization step 
%
% parameters are 
% p(1)=lqw,    log angular velocity variance
% p(2)=lr,     log measurement variance
% p(3:3+c-1)   log component variances
global dt

new = p;
for j=1:numel(gi)
    if(gi(j)==1) % dlb/dqw
        new(1) = (log(I2(1,1))-log(N))/2;
    end
    if(gi(j) >= 3) % dlb/dqx(ri)
      I = I2;
      c1 = 2*(I(2,2)+I(4,4)+I(6,6))/dt^3;
      c2 = I(2,3)+I(3,2)+I(4,5)+I(5,4)+I(6,7)+I(7,6);
      c2 = c2/dt^2;
      c3 = 2*(I(3,3)+I(5,5)+I(7,7))/(3*dt);
      new(3) = (log(c1-c2+c3)-log(N))/2;
      
    end
    if(gi(j)==2) % dlb/dr
        new(2) = (log(I3(1,1))-log(N))/2;
    end
end


end




