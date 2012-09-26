% Setup

global dt m0 P0 H A


N = 1500;
T = 14;
dt = T/N;


%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%
qx = 0.2;       % parameterization in standard deviations
qy = 0.1;
r =  log(0.1);


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
gis = [0 0 0 0 1 1;
       0 0 1 0 0 1;
       0 0 1 0 1 0;
       0 0 1 0 1 1;
      ];
logi = [1 1 1 1 0 0]; logi = logi > 0;

fn = '../data/Ballistic_%s%.0f_%.0f';
iters = [30 100;
         30 100;
         30 100;];
NNs = [30 30 30];
% iters = [10 10;
%          10  10;
%          10 10;];
% NNs = [1 1 1];

for i=1:size(gis,1)


gi = find(gis(i,:)>0);

min_iter_em =   iters(i,1);
max_iter_em =   iters(i,2);
min_iter_bfgs = iters(i,1);
max_iter_bfgs = iters(i,2);
NN = NNs(i);
est_em =   zeros(numel(gi),max_iter_em,NN);
lh_em =   zeros(max_iter_em,NN);
est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);
lh_bfgs =   zeros(max_iter_em,NN);

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
%   figure(1); clf;
%   plot(xs(1,:),xs(3,:),ys(1,:),ys(2,:),'kx');
%   axis equal; grid on;
%   pause
  
 

  % INITIAL POINT
  p0 = p_true;
  p0(gis(i,:)&~logi) = p0(gis(i,:)&~logi)*(4*rand-2); % 0.5-1.5 * true
  p0(gis(i,:)&logi) = p0(gis(i,:)&logi) + log(4*rand-2);
  
  % EM
  
  pt = p_true; pt(logi) = exp(pt(logi));
  cel = pNames(gi);cel(2,:)=num2cell(pt(gi));
  fprintf(1,'TRUE:    %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  p0t = p0; p0t(logi) = exp(p0t(logi));
  cel = pNames(gi);cel(2,:)=num2cell(p0t(gi));
  fprintf(1,'INITIAL: %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  tic;
  [opt,lh,vals] = Ballistic_EM(p0,gi,ys,[],[],max_iter_em,min_iter_em);
  tm = toc;
  
  optt = p0; optt(gi) = opt; optt(logi) = exp(optt(logi));
  cel = pNames(gi);cel(2,:)=num2cell(optt(gi)');
  fprintf(1,'EM %.0f: %s\n\n\n\n',k,sprintf('%s: %5.4f ',cel{:}));
  
  num = size(vals,2);
  est_em(:,1:num,k) = vals;
  lh_em(1:num,k) = lh;
  if num < max_iter_em
    est_em(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_em-num);
    lh_em(num+1:end,k) = repmat(lh(end),1,max_iter_em-num);
  end
  evals_em(1,k) = tm;
  
  % BFGS
  pt = p_true; pt(logi) = exp(pt(logi));
  cel = pNames(gi);cel(2,:)=num2cell(pt(gi));
  fprintf(1,'TRUE:    %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  p0t = p0; p0t(logi) = exp(p0t(logi));
  cel = pNames(gi);cel(2,:)=num2cell(p0t(gi));
  fprintf(1,'INITIAL: %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  tic;
  [opt,lh,vals,fcn_evals] = Ballistic_BFGS(p0,gi,ys,1e-60,1e-60,max_iter_bfgs,min_iter_bfgs);
  tm = toc;
  optt = p0; optt(gi) = opt; optt(logi) = exp(optt(logi));
  cel = pNames(gi);cel(2,:)=num2cell(optt(gi)');
  fprintf(1,'BFGS %.0f: %s\n\n\n\n',k,sprintf('%s: %5.4f ',cel{:}));
  num = size(vals,2);
  est_bfgs(:,1:num,k) = vals;
  lh_bfgs(1:num,k) = lh;
  if num < max_iter_bfgs
    est_bfgs(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_bfgs-num);
    lh_bfgs(num+1:end,k) = repmat(lh(end),1,max_iter_bfgs-num);
  end
  evals_bfgs(:,k) = [tm;fcn_evals];
end

%plot(max_iter_em,est_em')

cel = pNames(gi);
save(sprintf(fn,sprintf('%s_',pNames{gi}),NN,N));

end





