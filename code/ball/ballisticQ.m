function [ Q ] = ballisticQ( qx,qy )
global dt


%   +-                        -+
%   |       5       4       3  |
%   |  Qx dt   Qx dt   Qx dt   |
%   |  ------, ------, ------  |
%   |    20      8       6     |
%   |                          |
%   |       4       3       2  |
%   |  Qx dt   Qx dt   Qx dt   |
%   |  ------, ------, ------  |
%   |    8       3       2     |
%   |                          |
%   |       3       2          |
%   |  Qx dt   Qx dt           |
%   |  ------, ------,  Qx dt  |
%   |    6       2             |
%   +-                        -+

    Q1 = [dt^4/20  dt^3/8 dt^2/6;
          dt^3/8   dt^2/3 dt/2;
          dt^2/6   dt/2     1];
         
    Q = blkdiag(qx*dt*Q1,qy*dt*Q1);

end

