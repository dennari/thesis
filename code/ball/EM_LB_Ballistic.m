function [glb] = EM_LB_Ballistic(p,m0T,gi,N,I1,I2,I3)
% em_lb_harmonic - lower bound and its gradient 
%
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement variance
% gi        a vector of indices to p, specifying the ones that are varying

    global P0 g0x g0y dt
    
  Q = ballisticQ2D(p(3),p(4));
  R = p(5)*eye(size(I3,1));
  m0 = [0 p(1) 0 p(2)]';

% gradient

  glb = zeros(numel(gi));
  for j=1:numel(gi)
      % assume all zero
      dSig = zeros(size(P0)); 
      dQ = zeros(size(Q));dR = zeros(size(R));
      dmu = zeros(size(Q,1),1);

      if(gi(j)==1) % dlb/dqw
          dQ = zeros(size(Q));
          dQ(1,1) = 1;
      end
      if(gi(j)==3) % dlb/dqx
        q = p(3);
        I11 = I2(1,1);
        I12 = I2(1,2);
        I21 = I2(2,1);
        I22 = I2(2,2);
        
        t1 = log(2*I22*dt^2 - (I12 + I21)*3*dt + 6*I11);
        t2 = 3*log(dt)+2*log(q);%log(dt^3*q^2);
        glb(j) = exp(t1-t2)-N/q;
        %glb(j) = (3*(2*I11 - I21*dt))/(dt^3*q^2) - 1/q - (3*I12 - 2*I22*dt)/(dt^2*q^2);

        dQ = ballisticQ2D(1,0);
        
        
      end
      if(gi(j)==4) % dlb/dqy
          dQ = ballisticQ2D(0,1);
      end
      if(gi(j)==5) % dlb/dr
          dR = eye(size(R,1));
          r = p(5);
          glb(j) = I3(1,1)/(2*r^2) + I3(2,2)/(2*r^2) - N/r;
      end

      % x_0
      glb1 = P0\(dSig/P0*I1+2*dmu*(m0T-m0)'-dSig);
      % x_k|x_(k-1)
      glb2 = Q\(dQ/Q*I2-N*dQ);
      % y_k|x_k
      glb3 = R\(dR/R*I3-N*dR);

      %glb(j) = 0.5*(trace(glb1)+trace(glb2)+trace(glb3));
  end


end




