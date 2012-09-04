function x_ = sinusoid_A(x,dt)
% The signal part of A is
%   +-                         -+
%   |                sin(dt w)  |
%   |    cos(dt w),  ---------  |
%   |                    w      |
%   |                           |
%   |  -w sin(dt w), cos(dt w)  |
%   +-                         -+
	%A = cos(dt*w)*eye(2)+rot90(diag([-w*sin(dt*w) sin(dt*w)/w]),1);
 	%A = cos(dt*w)*eye(2)+rot90(diag([-sin(dt*w) sin(dt*w)]),1);
  %A = blkdiag(A,1);
  
  A = @(w) [(cos(dt*w)), (sin(dt*w)/w), 0;
          (-w*sin(dt*w)), (cos(dt*w)),  0;
          0,                0,          1];
	
  x_ = zeros(size(x));
  for k=1:size(x,2)
     x_(:,k) = A(max(abs(x(3,k)),1e-9))*x(:,k); 
  end
 
end