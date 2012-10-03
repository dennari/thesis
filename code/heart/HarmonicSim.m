
%% Setup
global dt c m0 P0 h f Jh Jf
N = 1500;
T = 15;
dt = T/N;
K = (0:N)*dt;

% the parameters of this model
lqx = log(1.2);    % log(sqrt) Dynamic model noise spectral density
lqw = log(0.8);   % log(sqrt) angular velocity noise variance
lr =  log(0.05);   % log(sqrt) measurement noise



c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x);
Jf = @(x) sinusoid_Jf(x);
Q = sinusoid_Q(lqw,lqx);              

if Q == diag(diag(Q))
  SQ = sqrt(Q);
else
  SQ = chol(Q,'lower');
end
R = sinusoid_R(lr);
SR = chol(R,'lower');
m0 = [1.8*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


X = zeros(xDim,N+1);
Y = zeros(1,N+1);

x = m0;
%x(1) = 1.8*2*pi;
X(:,1) = x;
Y(:,1) = h(x);

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  X(:,k) = x;
	Y(:,k) = mvnrnd(h(x),R)';
  
end

figure(1); clf;
plot(K,H*X,K,Y,'kx'); grid on;

%% Compute
% parameters are 
% p(1)=lqw,    log angular velocity std
% p(2)=lr,     log measurement std
% p(3)=lqx     log component std

pNames = {'lqw' 'lr' 'lqx'};
p_true = [lqw lr lqx];
gis = [1 0 1;
       1 1 1;];
    
    
fn = '../data/HarmonicSim_%s%.0f_%.0f';
iters = ones(2,2)*25;

NNs = [10 10];


for i=1:size(gis,1)


gi = find(gis(i,:)>0);

min_iter_em =   iters(i,1);
max_iter_em =   iters(i,2);
min_iter_bfgs = round(iters(i,1)/2);
max_iter_bfgs = round(iters(i,2)/2);
NN = NNs(i);
est_em =   zeros(numel(gi),max_iter_em,NN);
lh_em =   zeros(max_iter_em,NN);
times_em = lh_em;

est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);
lh_bfgs =   zeros(max_iter_em,NN);
times_bfgs = lh_bfgs;

evals_em = zeros(1,NN);
evals_bfgs = zeros(max_iter_em,NN);



for k=1:NN
 

  %%%%%%%%%%%%%%
  %% SIMULATE %%
  %%%%%%%%%%%%%%
  
  X = zeros(xDim,N+1);
  Y = zeros(1,N+1);

  x = m0;
  %x(1) = 0.5*2*pi;
  X(:,1) = x;
  Y(:,1) = h(x);

  for j=2:N+1
    x = mvnrnd(f(x),Q)';
    X(:,j) = x;
    Y(:,j) = mvnrnd(h(x),R)';
  end
  
  
  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)*(0.2+1.6*rand);
  
  %%%%%%%%%%%%
  % EM %%%%%%%
  %%%%%%%%%%%%
  
  % PRINT
  pt = exp(p_true);
  cel = pNames(gi);cel(2,:)=num2cell(pt(gi));
  fprintf(1,'TRUE:    %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  p0t = exp(p0);
  cel = pNames(gi);cel(2,:)=num2cell(p0t(gi));
  fprintf(1,'INITIAL: %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  tic;
  % RUN
  [opt,lh,vals,times] = Harmonic_EM(p0,gi,Y,[],[],max_iter_em,min_iter_em);
  tm = toc;
  % SAVE
  nn = numel(lh);
  est_em(:,1:nn,k) = vals;
  lh_em(1:nn,k) = lh;
  times_em(1:nn,k) = times;
  fprintf('EM time/iter: %.4f\n',sum(times)/sum(times>0));
  
  %num = size(vals,2);
  % PRINT
 
  %optt = p0; optt(gi) = opt; optt=exp(optt);
  %cel = pNames(gi);cel(2,:)=num2cell(optt(gi)');
  %fprintf(1,'EM %.0f: %s\n\n\n\n',k,sprintf('%s: %5.4f ',cel{:}));
  
  
  %%%%%%%%%%%%%%
  % BFGS%%%%%%%%
  %%%%%%%%%%%%%%
  
  %PRINT
  pt = exp(p_true);
  cel = pNames(gi);cel(2,:)=num2cell(pt(gi));
  fprintf(1,'TRUE:    %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  p0t = exp(p0);
  cel = pNames(gi);cel(2,:)=num2cell(p0t(gi));
  fprintf(1,'INITIAL: %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  % RUN
  [opt,lh,vals,funccount,times] = Harmonic_BFGS(p0,gi,Y,[],[],max_iter_bfgs,min_iter_bfgs);
  %tm = toc;
  fprintf('BFGS time/funcCount: %.4f\n',sum(times)/sum(funccount));

  % SAVE
  %num = size(vals,2);
  nn = numel(lh);
  est_bfgs(:,1:nn,k) = vals;
  lh_bfgs(1:nn,k) = lh;
  times_bfgs(1:nn,k) = times;

  evals_bfgs(1:nn,k) = funccount;
  %msg
  
  % PRINT
  optt = exp(p0); optt(gi) = exp(opt);
  cel = pNames(gi);cel(2,:)=num2cell(optt(gi)');
  fprintf(1,'BFGS %.0f: %s\n\n\n\n',k,sprintf('%s: %5.4f ',cel{:}));
  
  
end
%break
%plot(max_iter_em,est_em')

cel = pNames(gi);
save(sprintf(fn,sprintf('%s_',pNames{gi}),NN,N));

end



