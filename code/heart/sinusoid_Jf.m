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
global c  
% c is the number of components

   
  J = zeros(numel(x));
  J(1,1) = 1;
  for j=1:c
    jj = (j-1)*2;
    J(jj+2,[1 jj+2 jj+3]) = grad1(j,x(1),x(jj+2),x(jj+3));
    J(jj+3,[1 jj+2 jj+3]) = grad2(j,x(1),x(jj+2),x(jj+3));

  end
  
  
 
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