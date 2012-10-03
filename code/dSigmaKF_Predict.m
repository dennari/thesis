function [dm_,dP_] = dSigmaKF_Predict(m,m_,S,f,dm,dP,dQ,Jf,usig,w,du)
    if size(w,2) < 2
      w(:,2) = w(:,1);
    end
    wm = w(:,1);
    NS = numel(wm);
    
    if nargin < 11 || isempty(du)
      du = zeros(size(usig,1),1);
    end
    
    %P = S*S';
    %S = chol(P,'lower');
    
    
    [dS,~] = dchol(S*S',dP);
    
    dm_ = 0;
    Jff = zeros([size(usig,1) size(usig,1) size(usig,2)]);
    for j=1:NS
      sigj = m+S*usig(:,j);
      dsigj = dm+dS*usig(:,j);
      Jff(:,:,j) = Jf(sigj); 
      dm_ = dm_ + Jff(:,:,j)*dsigj;
    end
    dm_ = wm(1)*dm_+du; % weights assumed equal
    
    dP_ = 0;
    for j=1:NS
      sigj = m+S*usig(:,j);
      dsigj = dm+dS*usig(:,j);
      %ff = [sigj(1);
      %      Jff(2:end,2:end,j)*sigj(2:end,:)];
            
      dP1 = (Jff(:,:,j)*dsigj-dm_)*(f(sigj)-m_)';
      dP_ = dP_ + dP1 + dP1';
    end
    dP_ = wm(1)*dP_+dQ;
    
  
   
end




