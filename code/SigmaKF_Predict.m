function [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w)


%   if nargin < 5 || isempty(u)
%     u = zeros(size(m));
%   end
  if size(w,2) < 2
    w(:,2) = w(:,1);
  end
  wm = w(:,1);
  wp = w(:,2);
  
	NS = size(usig,2);
  sig=f(S*usig+repmat(m,1,NS));%+repmat(u,1,NS);
  m_ = sig*wm;
  
  [~,S_] = qr([(sig-repmat(m_,1,NS))*diag(wp) SQ]',0);
  S_ = S_'; % we want lower triangular

 
   
end




