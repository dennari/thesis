function [glb] = EM_LB_Ballistic(p,m0T,gi,N,I1,I2,I3,smk,smkk,Q,u)
% parameters are 
%
% p{1}=lqx,  x process std
% p{2}=lqy,  y process std
% p{3}=lr,   measurement std

 global dt A   


% gradient

  glb = zeros(numel(gi),1);
  for j=1:numel(gi)

      if(gi(j)==1) % dlb/dqx
        lq = p(1); 
        q = exp(2*lq);
        I11 = I2(1,1);
        I12 = I2(1,2);
        I21 = I2(2,1);
        I22 = I2(2,2);

        
        t1 = log(2*I22*dt^2 - (I12 + I21)*3*dt + 6*I11);
        t2 = 3*log(dt)+2*log(q);
        glb(j) = exp(t1-t2)*2*q-2*N;
        
      end
      if(gi(j)==2) % dlb/dqy
        lq = p(2);
        q = exp(2*lq);
        I11 = I2(1+2,1+2);
        I12 = I2(1+2,2+2);
        I21 = I2(2+2,1+2);
        I22 = I2(2+2,2+2);
        
        t1 = log(2*I22*dt^2 - (I12 + I21)*3*dt + 6*I11);
        t2 = 3*log(dt)+2*log(q);
        glb(j) = exp(t1-t2)*2*q-2*N;
        
      end
      if(gi(j)==3) % dlb/dr
          lr = p(3);
          r = exp(2*lr);
          glb(j) = (I3(1,1)/(2*r^2) + I3(2,2)/(2*r^2) - N/r)*r*2;
      end
      if(gi(j)==4) % dlb/dq (joint process noise)
        q = exp(2*p(1));

        I = I2;
        t1 = log(2*(I(2,2)+I(4,4))*dt^2 ...
             -3*(I(1,2) + I(2,1) + I(3,4) + I(4,3))*dt ...
             +6*(I(1,1)+I(3,3)));
        t2 = 3*log(dt)+2*log(q);
        glb(j) = exp(t1-t2)*2*q-4*N;
        
      end
      if(gi(j)==5 || gi(j)==6) % dlb/du
        ii = gi(j)-3+(gi(j)-5);
        e = zeros(4,1);
        e(ii) = dt;
        glb(j) = trace(Q\(e*smk'-e*smkk'*A'-N*e*u'));
        
      end
      
      


  end


end




