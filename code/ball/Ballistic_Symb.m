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

I = sym('I%d%d',[2 2]);
lbr = -0.5*trace(Rs\I)-0.5*N*log(det(Rs));
dlbr_r = simple(2*diff(lbr,r)*r);
solve(dlbr_r==0,r)


%% Q

I = sym('I%d%d',[4 4]);
syms qx qy positive
Qf = blkdiag(Qs/q*qx,Qs/q*qy);
lbq = -0.5*trace(I/Qf)-0.5*N*log(det(Qf));
%simple(diff(lbq,qx))
dlbq_qx = simple(2*diff(lbq,qx)*qx);
dlbq_qy = simple(2*diff(lbq,qy)*qy);

anss = solve(dlbq_qx==0,dlbq_qy==0,qx,qy)


%% u

I = sym('I%d%d',[4 4]);
syms qx qy positive
syms ux uy real
mk = sym('mk%d',[4 1]); sym(mk,'real');
mkk = sym('mkk%d',[4 1]); sym(mkk,'real');
A = blkdiag(As,As);

Qf = blkdiag(Qs/q*qx,Qs/q*qy);
u = [0 dt*ux 0 dt*uy]';
dux = diff(u,ux);
dQux = trace(Qf\(dux*mk'-dux*mkk'*A'-N*dux*u'));
duy = diff(u,uy);
dQuy = trace(Qf\(duy*mk'-duy*mkk'*A'-N*duy*u'));
%dlbq_q = simple(2*diff(lbq,qx)*qx);

anss = solve(dQux==0,dQuy==0,ux,uy)