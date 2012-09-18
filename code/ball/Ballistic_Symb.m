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
A = simple(expm(F*dt));
n   = size(F,1);
Phi = [F L*q*L'; zeros(n,n) -F'];
AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
Q   = simple(AB(1:n,:)/AB((n+1):(2*n),:));
R   = r*eye(2);

%% Compute the LB and M-steps in closed form

I = sym('I%d%d',[2 2]);
lbq = -0.5*trace(Q\I)-0.5*N*log(det(Q));
dlbq = diff(lbq,q)

lbr = -0.5*trace(R\I)-0.5*N*log(det(R));
dlbr = diff(lbr,r)

solve(dlbq==0,q)
solve(dlbr==0,r)