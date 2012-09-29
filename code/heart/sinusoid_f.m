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
global dt c
% c is the number of components


  A = @(w) [(cos(dt*w)), (sin(dt*w));
          -(sin(dt*w)), (cos(dt*w))];
   
    x_ = zeros(size(x));
    
    for k=1:size(x,2)
       w = x(1,k);
       x_(1,k) = w;%x(1,k); % noise driven
       for kk=1:c
         jj = (kk-1)*2; 
         x_(jj+2:jj+3,k) = A(kk*w)*x(jj+2:jj+3,k);
       end
    end
 
end