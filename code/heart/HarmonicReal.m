
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


ds = round(60*secs/60);
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
% break


% the parameters of this model
lqx = log(0.7);           % log Dynamic model noise spectral density
lqw = log(0.5);           % log angular velocity variance
lr =  log(0.0001);          % log measurement noise


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
gis = [1 0 1;
       1 1 1;];
    
    
logi = [1 1 1]; logi = logi > 0;

fn = '../data/Harmonic_%s%.0f_%.0f';
iters = ones(2,2)*50;

NNs = [5 5];


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
 

  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)*(0.9+0.2*rand); % 0.8-1.2 * true
  %p0(gis(i,:)&logi) = p0(gis(i,:)&logi) + log(4*rand-2);
  
  %%%%%%%%%%%%
  % EM %%%%%%%
  %%%%%%%%%%%%
  
  % PRINT
  pt = p_true; pt(logi) = exp(pt(logi));
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
  pt = p_true; pt(logi) = exp(pt(logi));
  cel = pNames(gi);cel(2,:)=num2cell(pt(gi));
  fprintf(1,'TRUE:    %s\n',sprintf('%s: %5.4f ',cel{:}));
  
  p0t = p0; p0t(logi) = exp(p0t(logi));
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
  optt = p0; optt(gi) = opt; optt(logi) = exp(optt(logi));
  cel = pNames(gi);cel(2,:)=num2cell(optt(gi)');
  fprintf(1,'BFGS %.0f: %s\n\n\n\n',k,sprintf('%s: %5.4f ',cel{:}));
  
  
end
%break
%plot(max_iter_em,est_em')

cel = pNames(gi);
save(sprintf(fn,sprintf('%s_',pNames{gi}),NN,N));

end


%% Plot DRIFTER frequency estimate

% S = load('../data/dataa_villelle.mat','card_freq','freq_t');
% 
% endi = find(S.freq_t<nm*60,1,'last');
% DRIFTER_f = S.card_freq(1:endi);
% DRIFTER_ft = S.freq_t(1:endi);
% plot(DRIFTER_ft,DRIFTER_f);




