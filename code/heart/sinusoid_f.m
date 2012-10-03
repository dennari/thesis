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
if nargin < 2
  global dt
end
% c is the number of components


  %A = @(w) [(cos(dt*w)), (sin(dt*w));
  %        -(sin(dt*w)), (cos(dt*w))];
        
A = @(x1)  [  cos(dt*x1), sin(dt*x1),             0,            0,             0,            0;
             -sin(dt*x1), cos(dt*x1),             0,            0,             0,            0;
                       0,          0,  cos(2*dt*x1), sin(2*dt*x1),             0,            0;
                       0,          0, -sin(2*dt*x1), cos(2*dt*x1),             0,            0;
                       0,          0,             0,            0,  cos(3*dt*x1), sin(3*dt*x1);
                       0,          0,             0,            0, -sin(3*dt*x1), cos(3*dt*x1)];
 

% A = @(x1)  [  cos(dt*x1), sin(dt*x1)/x1,             0,            0,             0,            0;
%              -x1*sin(dt*x1), cos(dt*x1),             0,            0,             0,            0;
%                        0,          0,  cos(2*dt*x1), sin(2*dt*x1)/(2*x1),         0,            0;
%                        0,          0, -2*x1*sin(2*dt*x1), cos(2*dt*x1),           0,            0;
%                        0,          0,             0,            0,  cos(3*dt*x1), sin(3*dt*x1)/(3*x1);
%                        0,          0,             0,            0, -3*x1*sin(3*dt*x1), cos(3*dt*x1)];
 
                           
                     
        
   
    x_ = x;
    for k=1:size(x,2)
      x_(2:end,k) = A(x(1,k))*x(2:end,k);
    end
 
end