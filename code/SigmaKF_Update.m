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

  [~,U] = qr([CW2 SR; CW1 zeros(XD,YD)]',0);
  U = U';

  Sy = U(1:YD,1:YD);
  %if Sy < 0
    %disp('J00');
    %disp(U);
  %end
  
  
  C = U((YD+1):end,1:YD);

  % the gain
  K = C/Sy;
  m = m+K*(y-my);
  S = U((YD+1):end,(YD+1):end);
  
  % propagate through the measurement function
  %sigp=h(sig);
  % apply the integration formula, wi(1,:) are the weights for mean
  %ym = sigp*w;
  %ym_rep = repmat(ym,1,NS);

  %d = sigp - ym_rep;
  %S = d*diag(w)*d'+SR;
  %C = (sig-m_rep)*diag(w)*(sigp-ym_rep)';

  % Kalman gain
  %K = C/S;
  % residual
  %d = y-ym;
  % updated mean
  %m = m+K*(y-my);
  % updated covariance matrix
  %S = S-K*S*K';
        
	





