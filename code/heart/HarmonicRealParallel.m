
%% Setup
global dt c m0 P0 h f Jh Jf

S = load('../data/dataa_villelle.mat','card_data','data_t');
y = S.card_data;
K = S.data_t;
clear S;
warning off;
offset = 20;
secs = 30;
%secs = 1;
starti = find(K<offset,1,'last'); 
endi = find(K<offset+secs,1,'last');
y = y(starti:endi);
K = K(starti:endi);
K = K-K(1);

ds = round(15*secs/60);
%ds = 15;
% downsample by ds
y = y(1:ds:end)';
K = K(1:ds:end);

dt = K(2)-K(1);
N = length(K);
T = K(end);

% plot first 60 seconds
% endi = find(K<5,1,'last');
% figure(1); clf;
% plot(K(1:endi),y(1:endi),'*-');
 %break


% the parameters of this model
lqx =-3.3; %log(0.2);           % log Dynamic model noise spectral density
lqw = -0.5;%log(0.2);           % log angular velocity variance
lr =  log(0.001);          % log measurement noise


c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x);
Jf = @(x) sinusoid_Jf(x);
m0 = [2*pi*1.1 zeros(1,xDim-1)]';
P0 = eye(xDim);


% adjust K and y to include the zeroth measurement
K = [0 K+dt];
Y = [H*m0 y];



% parameters are 
% p(1)=lqw,    log angular velocity std
% p(2)=lr,     log measurement std
% p(3)=lqx     log component std

pNames = {'lqw' 'lr' 'lqx'};
p_true = [lqw lr lqx];
gis = [1 0 1;]
%       1 1 1;];
    
    
fn = '../data/Harmonic_%s%.0f_%.0f';
iters = ones(2,2)*80;

NNs = [100 10];


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

jobs = cell(2,NN);

d = HarmonicDisp();
for k=1:NN
 

  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)*(0.1+2*rand); % 0.8-1.2 * true
  
  jobs{1,k} = ...
    batch(@myEM,4,{@Harmonic_LH,@EM_M_Harmonic,p0,gi,Y,@d.dispIter,[],[],...
    max_iter_em,min_iter_em,{dt c m0 P0 h f}},'PathDependencies','./heart');
 
  % RUN
  jobs{2,k} = ...
    batch(@myBFGS,5,{@Harmonic_LH,p0,gi,Y,@d.dispIter,[],[],...
    max_iter_bfgs,min_iter_bfgs,{dt c m0 P0 h f Jh Jf}},'PathDependencies','./heart');

  
end

end

%%
for k=1:NN
 
  
  %%%%%%%%%%%%
  % EM %%%%%%%
  %%%%%%%%%%%%
  

  % SAVE
  o = fetchOutputs(jobs{1,k});
  [~,lh,vals,times] = o{:};
  nn = numel(lh);
  est_em(:,1:nn,k) = vals;
  lh_em(1:nn,k) = lh;
  times_em(1:nn,k) = times;
  
  %%%%%%%%%%%%%%
  % BFGS%%%%%%%%
  %%%%%%%%%%%%%%
  

  % SAVE
  o = fetchOutputs(jobs{2,k});
  [~,lh,vals,funccount,times] = o{:};
  nn = numel(lh);
  est_bfgs(:,1:nn,k) = vals;
  lh_bfgs(1:nn,k) = lh;
  times_bfgs(1:nn,k) = times;
  evals_bfgs(1:nn,k) = funccount;
  
end


cel = pNames(gi);
save(sprintf(fn,sprintf('%s_',pNames{gi}),NN,N));

