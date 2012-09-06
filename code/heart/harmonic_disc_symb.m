%% Discretize stochastic oscillator in closed form

syms w dt g Qx real;

g = 0;
w = 0;
F = [ 0  1; 
    -w^2  -g];
L = [0;1];
A = simple(expm(F*dt));
n   = size(F,1);
Phi = [F L*Qx*L'; zeros(n,n) -F'];
AB  = expm(Phi*dt)*[zeros(n,n);eye(n)];
Q   = simple(AB(1:n,:)/AB((n+1):(2*n),:));

%% test

w_r = 2*pi*2;
dt_r = 0.02;
g_r = 0.4;
Qc_r = 0.01;

[Ar,Qr] = splti_disc(subs(F,[w g],[w_r g_r]),L,Qc_r,dt_r)
subs(A,[w g dt],[w_r g_r dt_r])
subs(Q,[g dt Qc],[g_r dt_r Qc_r])
