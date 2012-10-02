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
      new(1) = trans(I2(1,1)/(dt*N));
    end
    if(gi(j)== 3) % dlb/dqx
      c3 = (I2(3,3)+I2(5,5)+I2(7,7))/(3*dt*N);
      new(3) = trans(c3);
    end
    if(gi(j)==2) % dlb/dr
      new(2) = trans(I3(1,1)/N);
    end
end


end

function t=trans(x)
  t = log(x)/2;
end




