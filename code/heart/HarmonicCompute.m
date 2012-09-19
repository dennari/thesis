% Setup

global dt H c m0 P0 h f

S = load('../data/dataa_villelle.mat','card_data','data_t');
y = S.card_data;
K = S.data_t;
clear S;
nm = 1;
% take first nm minutes
endi = find(K<nm*60,1,'last');
y = y(1:endi);
K = K(1:endi);


ds = 30;
% downsample by ds
y = y(1:ds:end)';
K = K(1:ds:end);

dt = K(2)-K(1);
N = length(K);
T = K(end);


% plot first 60 seconds
%endi = find(K<60,1,'last');
%figure(1); clf;plot(K(1:endi),y(1:endi));
%break;

% the initial parameters
lqx = log(0.7);           % log Dynamic model noise spectral density
lqw = log(0.04);           % log angular velocity variance
lr =  log(0.02);          % log measurement noise



c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
f = @(x) sinusoid_f(x);
m0 = zeros(xDim,1);
m0(1) = 1.1*2*pi;
P0 = eye(xDim);

% adjust K and y to include the zeroth measurement
K = [0 K+dt];
Y = [H*m0 y];


% optimize


% parameters are 
% p(1)=lqw,    log angular velocity variance
% p(2)=lr,     log measurement variance
% p(3:3+c-1)   log signal component variances


p00 = [lqw lr repmat(lqx,1,c)];
gi = 1; % which one we're estimating


min_iter_em =   50;
max_iter_em =   50;
min_iter_bfgs = 10;
max_iter_bfgs = 10;
NN = 10;
est_em =   zeros(numel(gi),max_iter_em,NN);
est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);

evals_em = zeros(1,NN);
evals_bfgs = zeros(2,NN);

ys = zeros(2,N+1);

for k=1:NN
  
  % INITIAL POINT
  p0 = p00;
  p0(gi) = p0(gi)*(rand+0.5);
  
  % EM
  tic;
  [~,~,vals] = Harmonic_EM(p0,gi,Y,[],[],max_iter_em,min_iter_em);
  tm = toc;
  fprintf(1,'EM round %.0f time: %.2f s\n',k,tm);
  est_em(:,:,k) = vals;
  evals_em(1,k) = tm;
  
  % BFGS
  tic;
  [~,~,vals,fcn_evals] = Harmonic_BFGS(p0,gi,Y,[],[],max_iter_bfgs,min_iter_bfgs);
  tm = toc;
  fprintf(1,'BFGS round %.0f time: %.2f s\n',k,tm);
  num = size(vals,2);
  est_bfgs(:,1:num,k) = vals;
  if num < max_iter_bfgs
    est_bfgs(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_bfgs-num);
  end
  evals_bfgs(:,k) = [tm;fcn_evals];
end

%save(sprintf('../data/HarmonicResultsR%.0f.mat',NN),...
%     'est_em','evals_em','max_iter_em','p_true','est_bfgs','evals_bfgs','max_iter_bfgs');







