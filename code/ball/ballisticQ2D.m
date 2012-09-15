function [ Q ] = ballisticQ2D( qx,qy )
global dt

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

    Q1 = [dt^2/3  dt/2;
          dt/2    1];
         
    Q = blkdiag(qx*dt*Q1,qy*dt*Q1);

end

