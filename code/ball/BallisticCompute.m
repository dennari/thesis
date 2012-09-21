% Setup

global dt m0 P0 H A g0x g0y u


N = 700;
T = 14;
dt = T/N;


%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%
lqx = log(0.8);       % parameterization in logarithm of standard deviation
lqy = lqx-log(3);
lr =  log(3);


K = (0:N)*dt;
v0 = 300/3.6; % magnitude of the initial velocity


alpha0 = (60/180)*pi; % initial direction
v0x = v0*cos(alpha0);
v0y = v0*sin(alpha0);
g0y = -9.81; % initial y acceleration
g0x = -1; % initial x acceleration
A1 = [1 dt; 
      0 1];  
A = blkdiag(A1,A1); % for two dimensions
H = zeros(2,4); H(1,1) = 1; H(2,3) = 1;



Q = ballisticQ2D(lqx,lqy);
SQ = chol(Q,'lower');
h = @(x) H*x;
R = ballisticR(lr);
SR = chol(R,'lower');
u = [0 dt*g0x 0 dt*g0y]';


m0 = [0 v0x 0 v0y]';
P0 = eye(size(m0,1));%diag([1e-6 7^-2 1e-6 7^-2]);


% optimize
% 1=lqx,  x process log(std)
% 2=lqy,  y process log(std)
% 3=lr,   measurement log(std)
pNames = {'qx' 'qy' 'r'};
gis = eye(3);
fn = '../data/Ballistic_%s_%.0f';
iters = [100 100;
         10  10;
         100 100;];
NNs = [100 100 100];
% iters = [10 10;
%          10  10;
%          10 10;];
% NNs = [1 1 1];

for i=1:size(gis,1)


gi = find(gis(i,:)>0);
p_true = [lqx lqy lr]; % initial guess
true = p_true(gi);
min_iter_em =   iters(i,1);
max_iter_em =   iters(i,2);
min_iter_bfgs = iters(i,1);
max_iter_bfgs = iters(i,2);
NN = NNs(i);
est_em =   zeros(numel(gi),max_iter_em,NN);
est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);

evals_em = zeros(1,NN);
evals_bfgs = zeros(2,NN);

xs = zeros(4,N+1);
ys = zeros(2,N+1);

x0 = m0;%mvnrnd(m0,P0)';
xs(:,1) = x0;
ys(:,1) = H*x0;

fprintf(1,'ESTIMATING %.0f: %.3f\n',gi,exp(p_true(gi)));

for k=1:NN
  % SIMULATE
  x = x0;
  for j=2:N+1
    x = mvnrnd(A*x,Q)'+u;
    xs(:,j) = x;
    ys(:,j) = mvnrnd(H*x,R)';
  end
%   figure(1); clf;
%   plot(xs(1,:),xs(3,:),ys(1,:),ys(2,:),'kx');
%   axis equal; grid on;
%   pause
  
  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)+log(rand+0.5); % 0.5-1.5 * true
  
  % EM
  tic;
  [~,~,vals] = Ballistic_EM(p0,gi,ys,[],[],max_iter_em,min_iter_em);
  tm = toc;
  fprintf(1,'EM round %.0f time: %.2f s\n',k,tm);
  est_em(:,:,k) = vals;
  evals_em(1,k) = tm;
  
  % BFGS
  tic;
  [~,~,vals,fcn_evals] = Ballistic_BFGS(p0,gi,ys,[],[],max_iter_bfgs,min_iter_bfgs);
  tm = toc;
  fprintf(1,'BFGS round %.0f time: %.2f s\n',k,tm);
  num = size(vals,2);
  est_bfgs(:,1:num,k) = vals;
  if num < max_iter_bfgs
    est_bfgs(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_bfgs-num);
  end
  evals_bfgs(:,k) = [tm;fcn_evals];
end

%plot(max_iter_em,est_em')

save(sprintf(fn,pNames{i},NN));

end





