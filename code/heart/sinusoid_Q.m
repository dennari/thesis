function Q_ = sinusoid_Q(lqw,lqx)
global dt

% the parameterization is in logarithms of variance
qw = exp(lqw);
qx = exp(lqx);

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
     
    c = cell(1,numel(qx)+1);
    c{1} = qw;
    for k = 1:numel(qx)
        c{k+1} = qx(k)*dt*[Q11 Q12; Q12 Q22];
    end
    
    Q_ = blkdiag(c{:});
   
end

