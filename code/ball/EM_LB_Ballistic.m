function [glb] = EM_LB_Ballistic(p,m0T,gi,N,I1,I2,I3)
% em_lb_harmonic - lower bound and its gradient 
%
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement variance
% gi        a vector of indices to p, specifying the ones that are varying

    global P0 g0x g0y
    
  Q = ballisticQ(p(3),p(4));
  R = p(5)^2*eye(size(I3,1));
  m0 = [0 p(1) g0x 0 p(2) g0y]';

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
          dQ = ballisticQ(1,0);
      end
      if(gi(j)==4) % dlb/dqy
          dQ = ballisticQ(0,1);
      end
      if(gi(j)==5) % dlb/dr
          dR = eye(size(R,1));
      end

      % x_0
      glb1 = P0\(dSig/P0*I1+2*dmu*(m0T-m0)'-dSig);
      % x_k|x_(k-1)
      glb2 = Q\(dQ/Q*I2-N*dQ);
      % y_k|x_k
      glb3 = R\(dR/R*I3-N*dR);

      glb(j) = 0.5*(trace(glb1)+trace(glb2)+trace(glb3));
  end


end




