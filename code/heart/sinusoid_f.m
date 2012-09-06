function x_ = sinusoid_f(x,dt)
%sinusoid_f the f(x_(k-1)) function of the hoscillator SSM
% x - the evaluation points as [XD,n] matrix
% dt - the discretization delta t
% x_ - returns [XD,n] matrix

% The signal part of A is
%   +-                         -+
%   |                sin(dt w)  |
%   |    cos(dt w),  ---------  |
%   |                    w      |
%   |                           |
%   |  -w sin(dt w), cos(dt w)  |
%   +-                         -+
  
  
  A = @(w) [(cos(dt*w)), (sin(dt*w)/w);
          (-w*sin(dt*w)), (cos(dt*w))];
   
  c = (size(x,1)-1)/2; % number of harmonics    
  x_ = zeros(size(x));
  for k=1:size(x,2)
     %A_ = blkdiag(A(max(abs(x(3,k)),1e-9),c);
     As = cell(1,c);
     % frequency is the first component of the state
     w = abs(x(1,k));
     for j=1:c
         As{j} = A(j*w); % multiples of fundamental w
     end
     x_(1,k) = w; % noise driven
     x_(2:end,k) = blkdiag(As{:})*x(2:end,k); 
  end
 
end