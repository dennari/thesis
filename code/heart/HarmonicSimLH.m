
%% Setup
global dt c
N = 500;
T = 15;
dt = T/N;
K = (0:N)*dt;

% the parameters of this model
lqx = log(0.6);    % log(sqrt) Dynamic model noise spectral density
lqw = log(0.5);   % log(sqrt) angular velocity noise variance
lr =  log(0.0001);   % log(sqrt) measurement noise



c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x);
Jf = @(x) sinusoid_Jf(x);
Q = sinusoid_Q(lqw,lqx);              

if Q == diag(diag(Q))
  SQ = sqrt(Q);
else
  SQ = chol(Q,'lower');
end
R = sinusoid_R(lr);
SR = chol(R,'lower');
m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


X = zeros(xDim,N+1);
Y = zeros(1,N+1);

x = m0;
x(1) = 0.5*2*pi;
X(:,1) = x;
Y(:,1) = h(x);

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  %x(1) = fr(k);
  X(:,k) = x;
	Y(:,k) = mvnrnd(h(x),R)';
  
end

figure(1); clf;
plot(K,H*X,K,Y,'kx'); grid on;

%% Compute
Nqx = 20;
Nqw = 20;
alpha = 0.8;
qx_range = linspace(lqx+log(alpha),lqx-log(alpha),Nqx);
qw_range = linspace(lqw+log(alpha),lqw-log(alpha),Nqw);

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 2
end


[QX,QW] = meshgrid(qx_range,qw_range);



fn = '../data/HarmonicSimLH_%.0f_%.0f';

[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

Z = zeros(Nqw,Nqx);
%tic
for i=1:Nqx
  for j=1:Nqw   
    SQ = sqrt(sinusoid_Q(QW(j,i),QX(j,i),0,dt,c));
    lh = 0; m = m0; S = P0;
    %tmp = zeros(Nqw,1);
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);

        if k==N+1; break; end; 

      
        [m,S,K,my,CSy] = SigmaKF_Update(m_,S_,Y(:,k+1),h,SR,usig,w);
        %%% CSy and CC are Cholesky decompositions _HERE_ %%%
        Sy = CSy*CSy';
        lh = lh + likelihood(Y(:,k+1)-my,Sy);

        
    end
    %tmp(j) = lh;
    Z(j,i) = lh;
    fprintf('%.2f\n',100*((i-1)*Nqw+j)/(Nqx*Nqw));
  end
  
  
end
%tm = toc
save(sprintf(fn,Nqx,Nqw));

%%

surf(exp(QX),exp(QW),Z)



