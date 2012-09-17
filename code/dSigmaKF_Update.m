function [dm,dP,dmy,dSy] = dSigmaKF_Update(m_,S_,h,dm_,dP_,K,my,Sy,C,y,dR,Jh,usig,w)
    if size(w,2) < 2
      w(:,2) = w(:,1);
    end
    wm = w(:,1);
    NS = numel(wm);

    
    [dS_,~] = dchol(S_*S_',dP_);
    
    dmy = 0;
    for j=1:NS
      sigj = m_+S_*usig(:,j);
      dsigj = dm_+dS_*usig(:,j);
      dmy = dmy + Jh(sigj)*dsigj;
    end
    dmy = wm(1)*dmy; % weights assumed equal
    
    dSy = 0;
    for j=1:NS
      sigj = m_+S_*usig(:,j);
      dsigj = dm_+dS_*usig(:,j);
      dP1 = (Jh(sigj)*dsigj-dmy) * (h(sigj)-my)';
      dSy = dSy + dP1 + dP1';
    end
    dSy = wm(1)*dSy+dR;
    
    dC = 0;
    for j=1:NS
      sigj = m_+S_*usig(:,j);
      dsigj = dm_+dS_*usig(:,j);
      dC = dC + (dsigj-dm_)*(h(sigj)-my)'+(sigj-m_)*(Jh(sigj)*dsigj-dmy)';
    end
    dC = wm(1)*dC;
    
    dK = dC/Sy+C*(Sy\dSy/Sy);
    
    dm = dm_+dK*(y-my)-K*dmy;
    dP = dP_-dK*S*K'-K*dS*K'-K*S*dK';
  
   
end




