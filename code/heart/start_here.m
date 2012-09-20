
%% Setup
global dt H c m0 P0 h f
N = 1500;
T = 25;
dt = T/N;

K = (0:N)*dt;

% the parameters of this model
lqx = log(0.7);           % log Dynamic model noise spectral density
lqw = log(0.001);           % log angular velocity variance
lr =  log(0.02);          % log measurement noise


c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
f = @(x) sinusoid_f(x);
Q = sinusoid_Q(lqw,repmat(lqx,1,c));              
SQ = chol(Q,'lower');
R = sinusoid_R(lr);
SR = sqrt(R);
m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.5*ones(1,cp(1));
 L3 = 2.5*ones(1,N-cp(2)+1);
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
%u = zeros(size(x));
y = ys;
MM = zeros(xDim,N+1); MM(:,1) = m0;
SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
dm0 = zeros(xDim,1); dP0 = zeros(xDim);
m = m0; lh = 0;  S = P0;
dm = dm0; dP = dP0;

for k=1:(N+1)
    [m_,S_] = SigmaKF_Predict(m,S,f,SQ,[],usig,w);
    if k==N+1; break; end; 

    [m,S,~,IM,IS] = SigmaKF_Update(m_,S_,y(:,k+1),h,SR,usig,w);
    
    
    MM(:,k+1) = m;
    SS(:,:,k+1) = S;
    likelihood(y(:,k+1)-IM,IS)
    lh = lh + likelihood(y(:,k+1)-IM,IS);

end
[MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); 



figure(1); clf;
plot(K,H*xs,K,H*MM,K,ys,'kx'); grid on;
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
% p(1)=lqw,    log angular velocity variance
% p(2)=lr,     log measurement variance
% p(3:3+c-1)   log signal component variances


p0 = [lqw lr repmat(lqx,1,c)];
gi = 3; % which one we're estimating
true = exp(p0(gi));

NN = 25;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;




rnge = abs(true)-log(3);
as = linspace(true-1.8*rnge,true+0.2*rnge,NN);

p = p0;
for j=1:NN
    j
    p(gi) = as(j);
    
    [lh,glh,MM,SS] = Harmonic_LH(p,ys,gi);
    %lh = Harmonic_LH(p,ys);
    lhs(j) = lh;
    glhs(j) = glh;
    [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); % D = Smoother Gain
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,ys,MS,SM,DD);
    glbs(j) = EM_LB_Harmonic(p,MS(:,1),gi,N,I1,I2,I3);
end

n = 4; m= 1;
figure(1); clf;
subplot(n,m,1);
plot(exp(as),lhs'); grid on; title('likelihood'); hold on;
plot([true true],ylim,'-r');
subplot(n,m,2);
plot(exp(as),glhs); grid on; title('dSens'); hold on;
plot([true true],ylim,'-r');
subplot(n,m,3);
plot(exp(as),glbs); grid on; title('dEM'); hold on;
plot([true true],ylim,'-r');
subplot(n,m,4);
plot(exp(as(1:end-1)),diff(lhs)./diff(as)); grid on; title('dNum'); hold on;
plot([true true],ylim,'-r');

save('../data/simulateHeartR.mat','lhs','glhs','glbs');

%% Test EM and BFGS


% parameters are 
% p(1)=lqw,    log angular velocity variance
% p(2)=lr,     log measurement variance
% p(3:3+c-1)   log signal component variances


p00 = [lqw lr repmat(lqx,1,c)];
gi = 1; % which one we're estimating


min_iter_em =   10;
max_iter_em =   10;
min_iter_bfgs = 15;
max_iter_bfgs = 15;
NN = 1;
est_em =   zeros(numel(gi),max_iter_em,NN);
est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);

evals_em = zeros(1,NN);
evals_bfgs = zeros(2,NN);

Y = ys;

for k=1:NN
  k 
 
  
  % INITIAL POINT
  %p0 = p00;
  %p0(gi) = p0(gi)*(rand+0.5);
  
  % EM
%   tic;
%   [~,~,vals] = Harmonic_EM(p0,gi,Y,[],[],max_iter_em,min_iter_em);
%   tm = toc;
%   est_em(:,:,k) = vals;
%   evals_em(1,k) = tm;
  
  % BFGS
  tic;
  [~,~,vals,fcn_evals] = Harmonic_BFGS(p0,gi,Y,[],[],max_iter_bfgs,min_iter_bfgs);
  tm = toc;
  num = size(vals,2);
  est_bfgs(:,1:num,k) = vals;
  if num < max_iter_bfgs
    est_bfgs(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_bfgs-num);
  end
  evals_bfgs(:,k) = [tm;fcn_evals];
end


    
