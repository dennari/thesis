function [m,S,K,my,Sy,C] = SigmaKF_Update(m,S,y,h,SR,usig,w)

  % the number of sigma points
  NS = size(usig,2);
  XD = numel(m);
  YD = numel(y);
  if size(w,2) < 2
    w(:,2) = w(:,1);
  end
  wm = w(:,1);
  wp = w(:,2);
  
  % UPDATE, p(x_k|y_1:k)=N(m,P)

  %m_rep = repmat(m,1,NS); 
  %sig = m_rep+chol(S,'lower')*usig;
  
  sig = S*usig;
  % centered and weighted
  CW1 = sig*diag(wp);

  % propagate through the dynamics function
  sig=h(sig+repmat(m,1,NS));
  my = sig*wm;
  CW2 = (sig-repmat(my,1,NS))*diag(wp);

  U = qr_ckf([CW2 SR; CW1 zeros(XD,YD)]');

  Sy = U(1:YD,1:YD);

  C = U((YD+1):end,1:YD);

  % the gain
  K = C/Sy;
  m = m+K*(y-my);
  S = U((YD+1):end,(YD+1):end);
  
	





