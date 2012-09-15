% Setup

global dt m0 P0 H A g0x g0y u

N = 1500;
T = 14;
dt = T/N;



K = (0:N)*dt;
v0 = 300/3.6; % magnitude of the initial velocity
qx = 0.9^2;
qy = qx;%0;%qx/10;

alpha0 = (60/180)*pi; % initial direction
v0x = v0*cos(alpha0);
v0y = v0*sin(alpha0);
g0y = -9.81; % initial y acceleration
g0x = -1; % initial x acceleration

Q = ballisticQ2D(qx,qy);
A1 = [1 dt; 
      0 1];  
A = blkdiag(A1,A1); % for two dimensions
H = zeros(2,4); H(1,1) = 1; H(2,3) = 1;
r = (0.9)^2;
R = r*eye(2);
u = [0 dt*g0x 0 dt*g0y]';

m0 = [0 v0x 0 v0y]';
P0 = eye(size(m0,1));%diag([1e-6 7^-2 1e-6 7^-2]);


% optimize

gi = 5;
p_true = [v0x v0y qx qy r]; % initial guess
true = p_true(gi);
min_iter_em =   10;
max_iter_em =   10;
min_iter_bfgs = 15;
max_iter_bfgs = 15;
NN = 20;
est_em =   zeros(numel(gi),max_iter_em,NN);
est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);

evals_em = zeros(1,NN);
evals_bfgs = zeros(2,NN);

xs = zeros(4,N+1);
ys = zeros(2,N+1);

x0 = m0;%mvnrnd(m0,P0)';
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=1:NN
  
  % SIMULATE
  x = x0;
  for j=2:N+1
    x = mvnrnd(A*x,Q)'+u;
    xs(:,j) = x;
    ys(:,j) = mvnrnd(H*x,R)';
  end
  
  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)*(rand+0.5);
  
  % EM
  tic;
  [~,~,vals] = Ballistic_EM(p0,gi,ys,[],[],max_iter_em,min_iter_em);
  tm = toc;
  est_em(:,:,k) = vals;
  evals_em(1,k) = tm;
  
  % BFGS
  tic;
  [~,~,vals,fcn_evals] = Ballistic_BFGS(p0,gi,ys,[],[],max_iter_bfgs,min_iter_bfgs);
  tm = toc;
  num = size(vals,2);
  est_bfgs(:,1:num,k) = vals;
  if num < max_iter_bfgs
    est_bfgs(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_bfgs-num);
  end
  evals_bfgs(:,k) = [tm;fcn_evals];
end

save(sprintf('../data/BallisticResults%.0f.mat',NN),...
     'est_em','evals_em','est_bfgs','evals_bfgs');







