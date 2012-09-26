function [new] = EM_M_Ballistic(p,m0T,gi,N,I1,I2,I3,smk,smkk)
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
    

    if(gi(j)==1) % dlb/dqx        lq = p(1); 
%         %q = exp(2*p(1));
%         I11 = I2(1,1);
%         I12 = I2(1,2);
%         I21 = I2(2,1);
%         I22 = I2(2,2);
% 
%         
%         t1 = 2*I22*dt^2 - (I12 + I21)*3*dt + 6*I11;
%         %glb(j) = exp(t1-t2)*2*q-2*N;
%         new(1) = (log(t1)-3*log(dt)-log(N))/2;
          t1 = log(6*I2(1,1) - 3*I2(1,2)*dt - 3*I2(2,1)*dt + 2*I2(2,2)*dt^2);
          t2 = 3*log(dt)+log(N);
          new(1) = (t1-t2)/2;
        
    end
    if(gi(j)==2) % dlb/dqy
          I = I2(3:4,3:4);
          t1 = log(6*I(1,1) - 3*I(1,2)*dt - 3*I(2,1)*dt + 2*I(2,2)*dt^2);
          t2 = 3*log(dt)+log(N);
          new(2) = (t1-t2)/2;
    end
    if(gi(j)==3) % dlb/dr
      new(3) = (log((I3(1,1)+I3(2,2)))-log(2*N))/2;
    end
    if(gi(j)==5) % dlb/dux
      new(5) = (3*smkk(1) - 3*smk(1) + 2*dt*smk(2) + dt*smkk(2))/(2*N*dt^2);
    end
    if(gi(j)==6) % dlb/dux
      new(6) = (3*smkk(3) - 3*smk(3) + 2*dt*smk(4) + dt*smkk(4))/(2*N*dt^2);
    end
end
