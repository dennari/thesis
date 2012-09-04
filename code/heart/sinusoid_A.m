function A_ = sinusoid_A(w,dt)
% The signal part of A is
%   +-                         -+
%   |                sin(dt w)  |
%   |    cos(dt w),  ---------  |
%   |                    w      |
%   |                           |
%   |  -w sin(dt w), cos(dt w)  |
%   +-                         -+
	%A = cos(dt*w)*eye(2)+rot90(diag([-w*sin(dt*w) sin(dt*w)/w]),1);
 	A = cos(dt*w)*eye(2)+rot90(diag([-sin(dt*w) sin(dt*w)]),1);
  A_ = blkdiag(A,1);
 
end