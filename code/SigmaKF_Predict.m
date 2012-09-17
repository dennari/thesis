function [m_,S_,dm_,dP_] = SigmaKF_Predict(m,S,f,SQ,u,usig,w,dm,dP,Jf)

  if nargin < 10 || isempty(Jf)
    Jf = [];
  end
  if nargin < 9 || isempty(dP)
    dP = [];
  end
  if nargin < 8 || isempty(dm)
    dm = [];
  end
  if nargin < 5 || isempty(u)
    u = zeros(size(m));
  end
  if size(w,2) < 2
    w(:,2) = w(:,1);
  end
  wm = w(:,1);
  wp = w(:,2);
  
	NS = size(usig,2);
  sig=f(S*usig+repmat(m,1,NS))+repmat(u,1,NS);
  m_ = sig*wm;
  
  [~,S_] = qr([(sig-repmat(m_,1,NS))*diag(wp) SQ]',0);
  S_ = S_'; % we want lower triangular

  if ~isempty(dm) && ~isempty(dP) && ~isempty(Jf)
    
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
    dP_ = wm(1)*dP_;
    
  end
   
end




