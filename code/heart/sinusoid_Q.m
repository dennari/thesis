function Q_ = sinusoid_Q(qw,qx,notransform)
global dt
if nargin < 3
   notransform = 0;
end
qx = qx(1);
qw = qw(1);
% the parameterization is in log(sqrt(var))
if ~notransform % by default it is assumed that input is in log(std)
  qw = exp(2*qw);
  qx = exp(2*qx); % NOT CHANGED !!!!!
end

% The signal part of Q is    
%   +-                -+
%   |       3       2  |
%   |  Qx dt   Qx dt   |
%   |  ------, ------  |
%   |    3       2     |
%   |                  |
%   |       2          |
%   |  Qx dt           |
%   |  ------,  Qx dt  |
%   |    2             |
%   +-                -+ 
        
    %Q11 = (2*dt*w - sin(2*dt*w))/2;
    %Q12 = sin(dt*w)^2;
    %Q22 = (2*dt*w + sin(2*dt*w))/2;

    Q11 = dt^2/3;
    Q12 = dt/2;
    Q22 = 1;
     
 
    Q1 = qx*dt*[Q11 Q12; Q12 Q22];
    Q_ = blkdiag(qw,Q1,Q1);
   
end

