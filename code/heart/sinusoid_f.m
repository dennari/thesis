function x_ = sinusoid_f(x)
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
global dt
% c is the number of components


  %A = @(w) [(cos(dt*w)), (sin(dt*w));
  %        -(sin(dt*w)), (cos(dt*w))];
        
A = @(x1)  [  cos(dt*x1), sin(dt*x1),             0,            0,             0,            0;
             -sin(dt*x1), cos(dt*x1),             0,            0,             0,            0;
                       0,          0,  cos(2*dt*x1), sin(2*dt*x1),             0,            0;
                       0,          0, -sin(2*dt*x1), cos(2*dt*x1),             0,            0;
                       0,          0,             0,            0,  cos(3*dt*x1), sin(3*dt*x1);
                       0,          0,             0,            0, -sin(3*dt*x1), cos(3*dt*x1)];
 
      
        
   
    x_ = x;
    for k=1:size(x,2)
      x_(2:end,k) = A(x(1,k))*x(2:end,k);
    end
 
end