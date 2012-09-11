function [X,P,IM,IS,LH] = kfsr_update(X,P,y,H,R)

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end

  %
  % update step
  %
  IM = H*X;
  IS = (R*R' + H*(P*P')*H');
  
  K = (P*P')*H'/IS;
  X = X + K * (y-IM);
  yDim = size(y,1);
  [~,RR] = qr([R H*P; zeros(size(X,1),size(y,1)) P]');
  P = RR((yDim+1):end,(yDim+1):end)';
  
  if nargout > 4
    LH = gauss_pdf(y,IM,IS);
  end

