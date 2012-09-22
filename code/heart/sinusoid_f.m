function x_ = sinusoid_f(x,testw)
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
global dt c ww
% c is the number of components
if nargin < 2
  testw = 0;
end

  A = @(w) [(cos(dt*w)), (sin(dt*w)/w);
          (-w*sin(dt*w)), (cos(dt*w))];
   
    x_ = zeros(size(x));
    %   As = cell(1,c);
     % frequency is the first component of the state
    %w = x(1,1);
    %x_(1,:) = ww;
    %for j=1:c
    %     As{j} = A(j*ww); % multiples of fundamental w
    %end
    %AA = blkdiag(As{:});
    for k=1:size(x,2)
       %A_ = blkdiag(A(max(abs(x(3,k)),1e-9),c);
       if testw
         w = ww;
       else
         w = x(1,k);%ww;
       end
       
       x_(1,k) = w;%x(1,k); % noise driven
       for kk=1:c
         jj = (kk-1)*2; 
         x_(jj+2:jj+3,k) = A(kk*w)*x(jj+2:jj+3,k);
       end
    end
 
end