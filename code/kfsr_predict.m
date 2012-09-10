function [x,P] = kfsr_predict(x,P,A,Q,B,u)

  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    B = [];
  end
  if nargin < 6
    u = [];
  end
  
  %
  % Apply defaults
  %
  if isempty(A)
    A = eye(size(x,1));
  end
  if isempty(Q)
    Q = zeros(size(x,1));
  end
  if isempty(B) && ~isempty(u)
    B = eye(size(x,1),size(u,1));
  end

  %
  % Perform prediction
  %
  [~,RR] = qr([P'*A'; Q']);
  xDim = size(x,1);
  P = RR(1:xDim,1:xDim)';
  if isempty(u)
    x = A * x;
  else
    x = A * x + B * u;
  end

