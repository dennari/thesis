function [dm_,dP_] = dSigmaKF_Predict(m,m_,S,f,dm,dP,dQ,Jf,usig,w)
    if size(w,2) < 2
      w(:,2) = w(:,1);
    end
    wm = w(:,1);


    
    [dS,~] = dchol(S*S',dP);
    
    dm_ = 0;
    for j=1:NS
      sigj = m+S*usig(:,j);
      dsigj = dm+dS*usig(:,j);
      dm_ = dm_ + Jf(sigj)*dsigj;
    end
    dm_ = wm(1)*dm_; % weights assumed equal
    
    dP_ = 0;
    for j=1:NS
      sigj = m+S*usig(:,j);
      dsigj = dm+dS*usig(:,j);
      dP1 = (Jf(sigj)*dsigj-dm_) * (f(sigj)-m_)';
      dP_ = dP_ + dP1 + dP1';
    end
    dP_ = wm(1)*dP_+dQ;
    
  
   
end




