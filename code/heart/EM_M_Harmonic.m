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
    % assume all zero
    
    if(gi(j)==1) % dlb/dqw
        new(1) = log(I2(1,1))-log(N);
    end
    if(gi(j) >= 3) % dlb/dqx
        I = I2;
        if sum(gi>=3) > 1 % assume each component has independent variance
          jj = (gi(j)-3)*2;
          mat = [1 dt;dt dt^2].*[6 -3; -3 2].*I(jj+2:jj+3,jj+2:jj+3);
          new(gi(j)) = log(sum(sum(mat)))-3*log(dt)-log(N);%/(N*dt^3);
        else % assume components have shared variance
          ans2 = 0;
          for c=1:3
            jj = (c-1)*2;
            ans2 = ans2 + sum(sum([1 dt;dt dt^2].*[6 -3; -3 2].*I(jj+2:jj+3,jj+2:jj+3)));
          end

          new(gi(j)) = log(ans2)-3*log(dt)-log(3)-log(N);%(3*N*dt^3);
        end
        
    end
    if(gi(j)==2) % dlb/dr
        new(2) = log(I3)-log(N);
    end
end


end




