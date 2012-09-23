% Setup
global dt c m0 P0 h f Jh Jf
N = 800;
T = 25;
dt = T/N;
K = (0:N)*dt;

% the parameters of this model
lqx = log(0.6);    % log(sqrt) Dynamic model noise spectral density
lqw = log(0.04);   % log(sqrt) angular velocity noise variance
lr =  log(0.01);   % log(sqrt) measurement noise



c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x);
Jf = @(x) sinusoid_Jf(x);
Q = sinusoid_Q(lqw,lqx);              

SQ = chol(Q,'lower');
R = sinusoid_R(lr);
SR = chol(R,'lower');
m0 = [exp(2*lqw) zeros(1,xDim-1)]';
P0 = eye(xDim);

% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.8*ones(1,cp(1));
 L3 = 1.2*ones(1,N-cp(2)+1);
 x = K((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(K(cp(2))-K(cp(1))))*(x-K(cp(1)))+L1(1);
 fr = 2*pi*[L1 L2 L3];

m0(1) = fr(1); 
% parameters are 
% p(1)=lqw,    log sqrt angular velocity variance
% p(2)=lr,     log sqrt measurement variance
% p(3)=qx      log sqrt signal component variances

pNames = {'qw' 'r' 'qx'};
gis = [1 0 1;]


fn = '../data/Harmonic_%s_%.0f_%.0f';
iters = ones(3,2)*30;
NNs = [30 30 30];
% iters = [10 10;
%          10  10;
%          10 10;];
%NNs = [1 1 1];


for i=1:size(gis,1)


gi = find(gis(i,:)>0);
p_true = [lqw lr lqx]; % initial guess
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

xs = zeros(xDim,N+1);
ys = zeros(1,N+1);


x0 = m0;%mvnrnd(m0,P0)';
xs(:,1) = x0;
ys(:,1) = H*x0;
fprintf(1,'ESTIMATING %.0f: %.3f\n',gi,exp(p_true(gi)));

for k=1:NN
  % SIMULATE  
  x = x0;
  for j=2:N+1
    x = mvnrnd(f(x),Q)';
    x(1) = fr(j);
    xs(:,j) = x;
    ys(:,j) = mvnrnd(h(x),R)';
  end
%   figure(1); clf;
%   plot(K,H*xs,K,ys,'kx'); grid on;
%   pause
  
  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)+log(rand+0.5); % 0.5-1.5 * true
  
  % EM
  tic;
  [~,~,vals] = Harmonic_EM(p0,gi,ys,[],[],max_iter_em,min_iter_em);
  tm = toc;
  fprintf(1,'EM round %.0f time: %.2f s\n',k,tm);
  num = size(vals,2);
  est_em(:,1:num,k) = vals;
  if num < max_iter_em
    est_em(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_em-num);
  end
  evals_em(1,k) = tm;
  
  % BFGS
  tic;
  [~,~,vals,fcn_evals] = Harmonic_BFGS(p0,gi,ys,[],[],max_iter_bfgs,min_iter_bfgs);
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

save(sprintf(fn,pNames{gi},NN,N));

end





