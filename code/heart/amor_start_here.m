
%% Setup
global dt H c m0 P0 h f
N = 500;
T = 25;
dt = T/N;

K = (0:N)*dt;

% the parameters of this model
lqx = log(0.9);    % log(sqrt) Dynamic model noise spectral density
lqw = log(0.015);   % log(sqrt) angular velocity noise variance
lr =  log(0.05);   % log(sqrt) measurement noise


c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
f = @(x) sinusoid_f(x);
Q = sinusoid_Q(lqw,repmat(lqx,1,c));              
SQ = chol(Q,'lower');
R = sinusoid_R(lr);
SR = chol(R,'lower');
m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.2*ones(1,cp(1));
 L3 = 0.3*ones(1,N-cp(2)+1);
 x = K((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(K(cp(2))-K(cp(1))))*(x-K(cp(1)))+L1(1);
 fr = 2*pi*[L1 L2 L3];

%A = 0.5; 
%fr = A*sin(2*pi*0.1*K)+1.5*A;
%fr = ones(size(K))*2*pi*0.8;

[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
xs = zeros(xDim,N+1);
ys = zeros(1,N+1);

m0(1) = fr(1); 
x0 = m0;
x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  x(1) = fr(k);
  xs(:,k) = x;
	ys(:,k) = mvnrnd(h(x),R)';
  
end

%% Compute LH and gradient on grid

% parameters are 
% p(1)=lqw,    log sqrt angular velocity variance
% p(2)=lr,     log sqrt measurement variance
% p(3:3+c-1)   log sqrt signal component variances


p0 = [lqw lr repmat(lqx,1,c)];
gi = 1; % which one we're estimating
true = p0(gi);

NN = 125;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;




start = true + log(0.4);
endd = true - log(0.4);
as = linspace(start,endd,NN);
%as = log(linspace(0.06,0.09,NN));

p = p0;
for j=1:NN
    j
    p(gi) = as(j);
    
    [lh,glh,MM,SS,SQ] = Harmonic_LH(p,ys,gi);
    lhs(j) = lh;
    glhs(j) = glh;
    
    [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); % D = Smoother Gain
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,ys,MS,SM,DD);
    glbs(j) = EM_LB_Harmonic(p,MS(:,1),gi,N,I1,I2,I3);
end

save('../data/HarmonicTesting.mat');




    
