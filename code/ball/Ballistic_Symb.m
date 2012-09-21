%% Discretize ballistic in closed form, random walk acceleration

syms dt Qx real;

F = [ 0  1 0; 
      0  0 1;
      0  0 0];
L = [0;0;1];
A = simple(expm(F*dt));
n   = size(F,1);
Phi = [F L*Qx*L'; zeros(n,n) -F'];
AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
Q   = simple(AB(1:n,:)/AB((n+1):(2*n),:));

%% Discretize ballistic in closed form, white noise acceleration

syms dt q r N positive;

F = [ 0  1; 
      0  0];
L = [0;1];
As = simple(expm(F*dt));
n   = size(F,1);
Phi = [F L*q*L'; zeros(n,n) -F'];
AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
Qs   = simple(AB(1:n,:)/AB((n+1):(2*n),:));
Rs   = r*eye(2);

%% Compute the LB and M-steps in closed form


lbr = -0.5*trace(R\I)-0.5*N*log(det(R));
dlbr_r = simple(2*diff(lbr,r)*r);


%% Q

I = sym('I%d%d',[4 4]);
syms qx qy positive
Qf = blkdiag(Qs/q*qx,Qs/q*qx);
lbq = -0.5*trace(I/Qf)-0.5*N*log(det(Qf));
%simple(diff(lbq,qx))
dlbq_q = simple(2*diff(lbq,qx)*qx);


% solve(dlbq==0,q)
% solve(dlbr==0,r)

