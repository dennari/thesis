
%% Setup
global dt H c m0 P0
N = 500;
T = 25;
dt = T/N;

K = (0:N)*dt;

% the parameters of this model
qx = 0.7;              % Dynamic model noise spectral density
qw = 0.1;               % angular velocity variance
r = 0.02;                % measurement noise
R = r;
SR = sqrt(r);

c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];

h = @(x) H*x;
f = @(x) sinusoid_f(x);
Q = sinusoid_Q(qw,repmat(qx,1,c));              
SQ = chol(Q,'lower');

m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
xs = zeros(xDim,N+1);
ys = zeros(1,N+1);

m0(1) = fr(1); 
x0 = m0;
x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;


y = ys;
MM = zeros(xDim,N+1); MM(:,1) = m0;
SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
m = m0; lh = 0;  S = P0;
for k=1:(N+1)
    [m_,S_] = SigmaKF_Predict(m,S,f,SQ,[],usig,w);
    if k==N+1; break; end; 

    [m,S,~,IM,IS] = SigmaKF_Update(m_,S_,y(:,k+1),h,SR,usig,w);
    
    
    MM(:,k+1) = m;
    SS(:,:,k+1) = S;
    %likelihood(y(:,k+1)-IM,IS)
    lh = lh + likelihood(y(:,k+1)-IM,IS);

end
[MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); 



figure(1); clf;
plot(K,ys,K,H*MM,K,H*MS); grid on;
figure(2); clf;
m = 3;
subplot(m,1,1);
plot(K,sqrt(sum((MM-xs).^2)),K,sqrt(sum((MS-xs).^2))); grid on; title('Err');
subplot(m,1,2);
plot(K,xs(1,:),K,MM(1,:),K,MS(1,:)); grid on; title('Freq');
subplot(m,1,3);
plot(K,squeeze(abs(SS(1,1,:))),K,squeeze(abs(SM(1,1,:)))); grid on; title('Freq Std');



%% Compute LH and gradient on grid

% parameters are 
% p(1)=qw,    angular velocity variance
% p(2)=r,     measurement variance
% p(3:3+c-1)  the signal component variances


p0 = [qw r repmat(qx,1,c)];
p = p0;

NN = 50;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;


gi = 2;

as = linspace(0.5*p0(gi),1.5*p0(gi),NN);


for j=1:NN
    j
    p(gi) = as(j);
    
    [lh,glh,MM,SS] = Harmonic_LH(p,ys,gi);
    lhs(j) = lh;
    glhs(j) = glh;
    [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); % D = Smoother Gain
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,ys,MS,SM,DD);
    glbs(j) = EM_LB_Harmonic(p,MS(:,1),gi,N,I1,I2,I3);
end

n = 3; m= 1;
figure(1); clf;
subplot(n,m,1);
plot(as,lhs'); grid on;
subplot(n,m,2);
plot(as,glhs,as,glbs); grid on;
subplot(n,m,3);
plot(as,sqrt((glbs-glhs).^2)); grid on;

%save('../../data/simulateHeartR.mat','lhs','glhs','glbs');

    
