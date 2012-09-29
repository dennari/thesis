function Q_ = sinusoid_Q(qw,qx,notransform)
global dt c
qx = qx(1);
if nargin < 3
   notransform = 0;
end

% the parameterization is in log(sqrt(var))
if ~notransform % by default it is assumed that input is in log(std)
  qx = exp(2*qx);
  qw = exp(2*qw);
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

    Q11 = 0;
    Q12 = 0;
    Q22 = 1;

    Q1 = qx*dt*[Q11 Q12; Q12 Q22];
    Q_ = zeros(2*c+1);
    Q_(1,1) = qw;
    for k=1:c
      jj = (k-1)*2;
      Q_(jj+2:jj+3,jj+2:jj+3) = Q1;
    end
   
end

