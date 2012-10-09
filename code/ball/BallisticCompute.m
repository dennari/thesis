% Setup

global dt m0 P0 H A


%N = 500;
T = 10;
dt = 0.005;
N = round(T/dt);

%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%
qx = 0.4;       % std
qy = 0.1;      % std
r =  log(1.5);  % log(std)
g0y = -9.81; % initial y acceleration
g0x = -1.8; % initial x acceleration
v0 = 40; % magnitude of the initial velocity
alpha0 = (60/180)*pi; % initial direction


v0x = v0*cos(alpha0);
v0y = v0*sin(alpha0);

A1 = [1 dt; 
      0 1];  
A = blkdiag(A1,A1); % for two dimensions
H = zeros(2,4); H(1,1) = 1; H(2,3) = 1;
Q = ballisticQ2D(qx,qy);
SQ = chol(Q,'lower');
h = @(x) H*x;
R = ballisticR(r);
SR = chol(R,'lower');
u = ballisticU(g0x,g0y);
m0 = [0 v0x 0 v0y]';
P0 = eye(size(m0,1));%diag([1e-6 7^-2 1e-6 7^-2]);


% optimize
% 1=qx,  x process std
% 2=qy,  y process std
% 3=r,   measurement log(std)
% 4=qxy, shared process std
% 5=ux,    acceleration x-component
% 6=uy,    acceleration y-component 
pNames = {'qx' 'qy' 'r' 'qxy' 'ux' 'uy'};
p_true = [qx qy r 0 g0x g0y];
gis = [0 0 1 0 1 1];
    
    
logi = [1 1 1 1 0 0]; logi = logi > 0;

fn = '../data/Ballistic_%s%.0f_%.0f';
iters = [100 100];
NNs = 100;
% iters = [10 10;
%          10  10;
%          10 10;];
% NNs = [1 1 1];
%%
d = BallisticDisp();
for i=1:size(gis,1)


gi = find(gis(i,:)>0);

min_iter_em =   iters(i,1);
max_iter_em =   iters(i,2);
min_iter_bfgs = iters(i,1);
max_iter_bfgs = iters(i,2);
NN = NNs(i);
est_em =   zeros(numel(gi),max_iter_em,NN);
lh_em =   zeros(max_iter_em,NN);
times_em = lh_em;

est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);
lh_bfgs =   zeros(max_iter_em,NN);
times_bfgs = lh_bfgs;

evals_em = zeros(1,NN);
evals_bfgs = zeros(max_iter_em,NN);

xs = zeros(4,N+1);
ys = zeros(2,N+1);

x0 = m0;
ys(:,1) = H*x0;




for k=1:NN
  % SIMULATE
  x = x0;
  for j=2:N+1
    x = mvnrnd(A*x,Q)'+u;
    ys(:,j) = mvnrnd(H*x,R)';
    if x(3) < 0; break; end;
  end
  ys = ys(:,1:j);
  N = j-1;

  
 

  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)*(2*rand-1); % 0.5-1.5 * true
  
  %%%%%%%%%%%%
  % EM %%%%%%%
  %%%%%%%%%%%%
  
  % PRINT
  fprintf('\n\n\nEM\n\n\n');
  d.dispHeader(gi);
  % RUN
  [~,lh,vals,times] = ...
    myEM(@Ballistic_LH,@EM_M_Ballistic,p0,gi,ys,@d.dispIter,[],[],...
    max_iter_em,min_iter_em);
  d.dispTrue(p_true,gi);
  % SAVE
  nn = numel(lh);
  est_em(:,1:nn,k) = vals;
  lh_em(1:nn,k) = lh;
  times_em(1:nn,k) = times;
  
  %%%%%%%%%%%%%%
  % BFGS%%%%%%%%
  %%%%%%%%%%%%%%
  
  % PRINT
  fprintf('\n\n\nBFGS\n\n\n');
  d.dispHeader(gi);
  % RUN
  [~,lh,vals,funccount,times] = ...
    myBFGS(@Ballistic_LH,p0,gi,ys,@d.dispIter,1e-60,1e-60,...
    max_iter_bfgs,min_iter_bfgs);
  d.dispTrue(p_true,gi);
  % SAVE
  nn = numel(lh);
  est_bfgs(:,1:nn,k) = vals;
  lh_bfgs(1:nn,k) = lh;
  times_bfgs(1:nn,k) = times;
  evals_bfgs(1:nn,k) = funccount;
  

end


cel = pNames(gi);
save(sprintf(fn,sprintf('%s_',pNames{gi}),NN,N));

end





