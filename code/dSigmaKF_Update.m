function [dm,dP,dmy,dSy] = dSigmaKF_Update(m_,S_,h,dm_,dP_,K,my,Sy,y,dR,Jh,usig,w)
    if size(w,2) < 2
      w(:,2) = w(:,1);
    end
    wm = w(:,1);
    NS = numel(wm);

    %P_ = S_*S_';
    %S_ = chol(P_,'lower');
    
    [dS_,~] = dchol(S_*S_',dP_);
    
    dmy = 0;
    Jhh = zeros([size(my,1) NS]);
    for j=1:NS
      sigj = m_+S_*usig(:,j);
      dsigj = dm_+dS_*usig(:,j);
      Jhh(:,j) = Jh(sigj)*dsigj;
      dmy = dmy + Jhh(:,j);
    end
    dmy = wm(1)*dmy; % weights assumed equal
    
    dSy = 0;
    dC = 0;
    for j=1:NS
      sigj = m_+S_*usig(:,j);
      dsigj = dm_+dS_*usig(:,j);
      d1 = (Jhh(:,j)-dmy);
      d2 = (h(sigj)-my);
      dP1 = d1 * d2';
      dSy = dSy + dP1 + dP1';
      dC = dC + (dsigj-dm_)*d2'+(sigj-m_)*d1';
      
    end
    dSy = wm(1)*dSy+dR;
    dC = wm(1)*dC;
    
    
    
    
    % this should not need to be computed again
    C = (S_*usig)*diag(wm)*(h(repmat(m_,1,NS)+S_*usig)-repmat(my,1,NS))';
    
    dK = dC/Sy-C*(Sy\dSy/Sy);
    
    dm = dm_+dK*(y-my)-K*dmy;
    dP = dP_-dK*Sy*K'-K*dSy*K'-K*Sy*dK';
  
   
end




