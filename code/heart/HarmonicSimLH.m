
%% Setup
%global dt c
N = 1500;
T = 15;
dt = T/N;
K = (0:N)*dt;

% the parameters of this model
qx = 0.4;    % log(sqrt) Dynamic model noise spectral density
qw = 0.2;   % log(sqrt) angular velocity noise variance
lr =  log(0.001);   % log(sqrt) measurement noise



c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x,dt);
Jf = @(x) sinusoid_Jf(x);

Qf = @(qw,qx) diag([qw^2, 0, qx^2, 0, qx^2, 0, qx^2]);

Q = Qf(qw,qx);
SQ = sqrt(Q);
R = sinusoid_R(lr);
SR = chol(R,'lower');
m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


X = zeros(xDim,N+1);
Y = zeros(1,N+1);

x = m0;
X(:,1) = x;
Y(:,1) = h(x);

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  %x(1) = fr(k);
  X(:,k) = x;
	Y(:,k) = mvnrnd(h(x),R)';
  
end

%figure(1); clf;
%plot(K,H*X,K,Y,'kx'); grid on;

%% Compute
Nqx = 100;
Nqw = 100;
alpha = 0.4;
qx_range = linspace(qx*alpha,qx/alpha,Nqx);
qw_range = linspace(qw*alpha,qw/alpha,Nqw);

if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 30
end


[QX,QW] = meshgrid(qx_range,qw_range);



fn = '../data/HarmonicSimLH_%.0f_%.0f';

[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

Z = zeros(Nqx,Nqw);
%tic
parfor i=1:Nqx
  tmp = zeros(Nqw,1);
  for j=1:Nqw   
    SQ = sqrt(Qf(QW(j,i),QX(j,i)));

    
    lh = 0; m = m0; S = P0;
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);

        if k==N+1; break; end; 
        yy = Y(k+1);
      
        [m,S,K,my,CSy] = SigmaKF_Update(m_,S_,yy,h,SR,usig,w);
        %%% CSy and CC are Cholesky decompositions _HERE_ %%%
        Sy = CSy*CSy';
        lh = lh + likelihood(yy-my,Sy);

        
    end
    tmp(j) = lh;
    
    %fprintf('%.2f\n',100*((i-1)*Nqw+j)/(Nqx*Nqw));
  end
  i
  Z(:,i) = tmp;
  
end
%tm = toc
save(sprintf(fn,Nqx,Nqw));

%%
load('tritonwrk/HarmonicSimLH_100_100.mat');
surf(QX,QW,Z); hold on;
plot3([qx qx],[qw qw],[min(min(Z)) max(max(Z))]);
xlabel('qx');
ylabel('qw');


