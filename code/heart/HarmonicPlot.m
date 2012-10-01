%% r - estimates
funs = plotFuns();
load('../data/Harmonic_lqw_lr_lqx_5_1501.mat');
NN = 5;
itr = max_iter_em;
times_em = times_em(1:itr,1:NN);
zerse = times_em <= 0;
% start from zero
times_em = times_em(1:itr,1:NN)-repmat(times_em(1,1:NN),itr,1);
lh_em = lh_em(1:itr,1:NN);
evals_em = cumsum(ones(itr,NN));
evals_em(lh_em <= 0) = 0;
est_em = est_em(:,1:itr,1:NN);
[lh_em_n,est_em_n] = funs.normalizeBFGS(evals_em,lh_em,est_em);



itr = max_iter_bfgs;
times_bfgs = times_bfgs(1:itr,1:NN);
zers = times_bfgs <= 0;
times_bfgs = times_bfgs-repmat(times_bfgs(1,1:NN),itr,1);
lh_bfgs = lh_bfgs(1:itr,1:NN);
est_bfgs = est_bfgs(:,1:itr,1:NN);
evals_bfgs = evals_bfgs(1:itr,1:NN);

[lh_bfgs_n,est_bfgs_n] = funs.normalizeBFGS(evals_bfgs,lh_bfgs,est_bfgs);
% LH
figure(1); clf;
% normalized
subplot(2,1,1);
y = lh_em_n;
plot(y);
%xlim([0 30]); ylim([18000 24500]);

subplot(2,1,2);
y = lh_bfgs_n;
plot(y);
%break;

%xlim([0 50]); ylim([-1e5 8000]);

% original
%subplot(2,1,2);
%times_bfgs(zers) = nan; lh_bfgs(zers) = nan;
%times_em(zerse) = nan; lh_em(zerse) = nan;
%plot(times_em,lh_em,'-b',times_bfgs,lh_bfgs,'-r');
%plot(times_bfgs/3600,lh_bfgs,'-r');

%xlim([0 3000]); ylim([22000 24500]);
%break
% mean over runs
%subplot(3,1,3);
%plot((0:itr-1)*avg_time_em,mean(lh_em,2),'-b',(0:itr-1)*avg_time_bfgs,mean(lh_bfgs_n,2),'-r');
%xlim([0 50]); ylim([-1e5 8000]);

% EST
figure(2); clf;
subplot(3,1,1);
x = 1:max_iter_em;
y = squeeze(est_em_n(1,x,:));
plot(x,y,'-b'); grid on;

subplot(3,1,2);
y = squeeze(est_em_n(2,x,:));
plot(x,y,'-b'); grid on;

subplot(3,1,3);
y = squeeze(est_em_n(3,x,:));
plot(x,y,'-b'); grid on;

figure(3); clf;
% original
subplot(3,1,1);
x = 1:max_iter_bfgs;
y = squeeze(est_bfgs_n(1,x,:));
plot(x,y,'-r'); grid on;

subplot(3,1,2);
y = squeeze(est_bfgs_n(2,x,:));
plot(x,y,'-r'); grid on;

subplot(3,1,3);
y = squeeze(est_bfgs_n(3,x,:));
plot(x,y,'-r'); grid on;

% figure(4); clf;
% % original
% subplot(3,1,1);
% x = 3:30;
% y1 = exp(squeeze(est_em(1,x,:))/log(10));
% y2 = exp(squeeze(est_bfgs_n(1,x,:))/log(10));
% semilogy(x,y1,'b',x,y2,'-r'); grid off;
% 
% subplot(3,1,2);
% y1 = exp(squeeze(est_em(2,x,:))/log(10));
% y2 = exp(squeeze(est_bfgs_n(2,x,:))/log(10));
% semilogy(x,y1,'b',x,y2,'-r'); grid off;
% 
% subplot(3,1,3);
% y1 = exp(squeeze(est_em(3,x,:))/log(10));
% y2 = exp(squeeze(est_bfgs_n(3,x,:))/log(10));
% semilogy(x,y1,'b',x,y2,'-r'); grid off;

% mean over runs
%subplot(2,1,2);
%plot((0:itr-1)*avg_time_em,mean(est1,2),'-b',(0:itr-1)*avg_time_bfgs,mean(est2n,2),'-r');
%xlim([0 50])

%% Export

%% trajectory

% figure(1); clf;
endi = find(K<7,1,'last');
% plot(K(1:endi),Y(1:endi),'k-x');
% ylim([-0.5 0.5]);
% break;

plt = struct();kw=struct(); kw1 = struct();
kw.alpha = 1.0; kw.ms = 4; kw.mfc = 'black';
kw1.lw = 0.9; kw1.alpha = 0.8;
plt.data = {{K(1:endi) Y(1:endi) '' kw1},{{K(1:endi) Y(1:endi) '*' kw}}};
plt.xlabel = '$t$';
%plt.legend = {'ECG' ''};
%plt.legendkw = struct('loc','lower center');
plt.w = textwidth*0.7;
plt.margins = [0.0 0.1 0.4 0.35];

%plotstruct(ax,plt);
pyplot('../img/harmonic_trajectory.pdf',plt);


%% Export

textwidth = 426.79134/72.27; % latex textwidth in inches
% plot the true locations and the measurements

% Likelihood
plt = struct();kw1=struct();
kw1.color = '#348ABD'; kw1.alpha=0.8; kw1.lw = 1.2;
x = 1:max_iter_em;
y = lh_em_n;
plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5+0.4;%/1.8;
plt.ticklabels = [0 1];
plt.margins = [0.0 0.0 0.1 0.5];
%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
yl = [1500 4500];
plt.axis = [min(x) max(x) yl];
pyplot('../img/harmonic_em_lh.pdf',plt,'../img/harmonic_em_lh.mat');

plt = struct();kw1=struct();
kw1.color = '#E24A33'; kw1.alpha=0.8; kw1.lw = 1.2;
x = 1:max_iter_bfgs;
y = lh_bfgs_n;
plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.1];
%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
plt.axis = [min(x) max(x) yl];
pyplot('../img/harmonic_bf_lh.pdf',plt,'../img/harmonic_bf_lh.mat');


%% Estimates
plt = struct();kw1=struct();
kw1.color = '#348ABD'; kw1.alpha=0.9; kw1.lw = 1.0;
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5+0.4;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.50];
x = 1:max_iter_em;
labelNames = {'\log(\sigma_\omega)', '\log(\sigma_r)','\log(\sigma_x)'};
yl = [min(x) max(x) -4 2;
      min(x) max(x) -25 -5;
      min(x) max(x) 0 -20];
for k = 1:3
  y1 = squeeze(est_em(k,:,:));
  plt.data = {{x y1 '' kw1},{' '}};
  var = pNames{gi(k)}; 
  plt.ylabel = sprintf('$%s$',labelNames{gi(k)});
  plt.axis = yl(gi(k),:);
  if gi(k) == 2 % r is last
    plt.margins(3) = 0.2;
    plt.ticklabels = [1 1];
  else
    plt.margins(3) = 0.1;
    plt.ticklabels = [0 1];
  end
  pyplot(sprintf('../img/harmonic_em_%s.pdf',var),plt,...
         sprintf('../img/harmonic_em_%s.mat',var));
end

plt = struct();kw1=struct();
kw1.color = '#E24A33'; kw1.alpha=0.8; kw1.lw = 1.2;
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.1];
x = 1:max_iter_bfgs;
yl = [min(x) max(x) -4 2;
      min(x) max(x) -25 -5;
      min(x) max(x) 0 -20];
for k = 1:3
  y1 = squeeze(est_bfgs_n(k,:,:));
  plt.data = {{x y1 '' kw1},{' '}};
  var = pNames{gi(k)};
  plt.axis = yl(gi(k),:);
  if gi(k) == 2 % r is last
    plt.margins(3) = 0.2;
    plt.ticklabels = [1 0];
  else
    plt.margins(3) = 0.1;
    plt.ticklabels = [0 0];
  end
  %plt.ylabel = sprintf('$%s$',labelNames{k});
  pyplot(sprintf('../img/harmonic_bf_%s.pdf',var),plt,...
         sprintf('../img/harmonic_bf_%s.mat',var));
end




