function [glb] = EM_LB_Ballistic(p,m0T,gi,N,I1,I2,I3)
% parameters are 
%
% p{1}=lqx,  x process std
% p{2}=lqy,  y process std
% p{3}=lr,   measurement std

 global dt
    


% gradient

  glb = zeros(numel(gi));
  for j=1:numel(gi)

      if(gi(j)==1) % dlb/dqx
        lq = p(1); 
        q = exp(2*lq);
        I11 = I2(1,1);
        I12 = I2(1,2);
        I21 = I2(2,1);
        I22 = I2(2,2);
        
        %t1 = log(2*I22*dt^2 - (I12 + I21)*3*dt + 6*I11);
        %t2 = 3*log(dt)+2*log(q);%log(dt^3*q^2);
        %glb(j) = exp(t1-t2)-N/q*exp(-p(gi(j)));
        glb(1) = (4*I22*dt^2 - 6*(I12 + I21)*dt + 12*I11)/(dt^3*q) - 2*N;
        
      end
      if(gi(j)==2) % dlb/dqy
        lq = p(2);
        q = exp(2*lq);
        I11 = I2(1+2,1+2);
        I12 = I2(1+2,2+2);
        I21 = I2(2+2,1+2);
        I22 = I2(2+2,2+2);
        
        %t1 = log(2*I22*dt^2 - (I12 + I21)*3*dt + 6*I11);
        %t2 = 3*log(dt)+2*log(q);%log(dt^3*q^2);
        %glb(j) = exp(t1-t2)-N/q*exp(-p(gi(j)));
        glb(2) = (4*I22*dt^2 - 6*(I12 + I21)*dt + 12*I11)/(dt^3*q) - 2*N;
      end
      if(gi(j)==3) % dlb/dr
          lr = p(3);
          r = exp(2*lr);
          glb(j) = (I3(1,1)/(2*r^2) + I3(2,2)/(2*r^2) - N/r)*r*2;
      end


  end


end




