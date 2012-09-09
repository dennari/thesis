%% Discretize ballistic in closed form

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

