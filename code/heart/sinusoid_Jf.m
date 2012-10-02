function J = sinusoid_Jf(x)
%sinusoid_Jf, Jacobian matrix of f
% x - the evaluation point as [XD,1] vector
% J - returns [XD,XD] matrix

% The signal part of A is
%   +-                         -+
%   |                sin(dt w)  |
%   |    cos(dt w),  ---------  |
%   |                    w      |
%   |                           |
%   |  -w sin(dt w), cos(dt w)  |
%   +-                         -+
global dt  c
% c is the number of components

   
  J = zeros(numel(x));
  J(1,1) = 1;
  for j=1:c
    jj = (j-1)*2;
    J(jj+2,[1 jj+2 jj+3]) = grad1(j,x(1),x(jj+2),x(jj+3));
    J(jj+3,[1 jj+2 jj+3]) = grad2(j,x(1),x(jj+2),x(jj+3));

  end
% x1 = x(1); x2=x(2); x3=x(3); x4= x(4); x5=x(5); x6 = x(6); x7= x(7);  
%   
%   J =[                                       1,           0,          0,             0,            0,             0,            0;
%            dt*x3*cos(dt*x1) - dt*x2*sin(dt*x1),  cos(dt*x1), sin(dt*x1),             0,            0,             0,            0;
%          - dt*x2*cos(dt*x1) - dt*x3*sin(dt*x1), -sin(dt*x1), cos(dt*x1),             0,            0,             0,            0;
%    2*dt*x5*cos(2*dt*x1) - 2*dt*x4*sin(2*dt*x1),           0,          0,  cos(2*dt*x1), sin(2*dt*x1),             0,            0;
%  - 2*dt*x4*cos(2*dt*x1) - 2*dt*x5*sin(2*dt*x1),           0,          0, -sin(2*dt*x1), cos(2*dt*x1),             0,            0;
%    3*dt*x7*cos(3*dt*x1) - 3*dt*x6*sin(3*dt*x1),           0,          0,             0,            0,  cos(3*dt*x1), sin(3*dt*x1);
%  - 3*dt*x6*cos(3*dt*x1) - 3*dt*x7*sin(3*dt*x1),           0,          0,             0,            0, -sin(3*dt*x1), cos(3*dt*x1)];
%   
%   
  
 
end

function g=grad1(j,x1,x2,x3)
  global dt
  g(1) = j*dt*x3*cos(dt*j*x1) - j*dt*x2*sin(dt*j*x1);
  g(2) = cos(j*dt*x1);
  g(3) = sin(j*dt*x1);
end
function g=grad2(j,x1,x2,x3)
  global dt;
  g(1) = -j*dt*x2*cos(dt*j*x1) - dt*j*x3*sin(dt*j*x1);
  g(2) = -sin(dt*j*x1);
  g(3) = cos(dt*j*x1);
end   