%% r - estimates
funs = plotFuns();
load('../data/Harmonic_lqw_lqx_100_3751.mat');
%NN = 5;
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
plot(y,'b-');
%xlim([0 30]); ylim([18000 24500]);

subplot(2,1,2);
y = lh_bfgs_n;
plot(y,'r');
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
for k=1:numel(gi)
  subplot(numel(gi),1,k);
  x = 1:max_iter_em;
  y = squeeze(est_em_n(k,x,:));
  plot(x,y,'-b'); grid on;
  title(pNames{gi(k)});
end

figure(3); clf;
for k=1:numel(gi)
  subplot(numel(gi),1,k);
  x = 1:max_iter_bfgs;
  y = squeeze(est_bfgs_n(k,x,:));
  plot(x,y,'-r'); grid on;
  title(pNames{gi(k)});
end

%% Plot3
figure(1); clf;
load('../data/HarmonicLH_60_60.mat');
mesh(QX,QW,Z); hold on;
% QX(1,48)
% QW(29,1)
Z1=Z;
QX1 = QX;
QW1 = QW;

%figure(2); clf;
load('../data/HarmonicLB_171_95.mat');
% QW(45,1)
% QX(1,58)
mesh(QX,QW,Z-(Z(50,50)-Z1(29,48))-0.03e4); hold on;

load('../data/HarmonicLB_116_95.mat');
% QW(45,1)
% QX(1,58)
mesh(QX,QW,Z-(Z(50,50)-Z1(29,48))+0.06e4); hold on;

xlim([-3.9  -1.5]);
ylim([-2.25 -0.1]);
zlim([1.6e4 2.0e4]);
break
%%
figure(1); clf;
plot3(squeeze(est_em_n(2,:,:)),squeeze(est_em_n(1,:,:)),...
      lh_em_n(:,:),'b*-'); grid on;
xlabel(pNames{gi(2)});ylabel(pNames{gi(1)});zlabel('lh');
zlim([1.6e4 2.0e4]); axis vis3d;

figure(2); clf;
plot3(squeeze(est_bfgs_n(2,:,:)),squeeze(est_bfgs_n(1,:,:)),...
      lh_bfgs_n(:,:),'r*-'); grid on;    
xlabel(pNames{gi(2)});ylabel(pNames{gi(1)});zlabel('lh');
zlim([1.6e4 2.0e4]); axis vis3d;



% est_em_n(:,10,6)

%figure(2); clf;
%plot3(squeeze(est_bfgs_n(1,:,:)),squeeze(est_bfgs_n(2,:,:)),...
%      lh_bfgs_n(:,:),'r*-'); grid on;
%xlabel(pNames{gi(1)});ylabel(pNames{gi(2)}); zlabel('lh');
%zlim([1.85e4 2.0e4]);

save('../img/harmonic_3D.mat','est_em_n','lh_em_n','est_bfgs_n','lh_bfgs_n','-v7');

%%


figure(1); clf;
subplot(1,2,1);
plot(squeeze(est_em_n(1,:,:)),...
      lh_em_n(:,:),'b*-'); grid on;
xlabel(pNames{gi(1)});ylabel('lh');
xlim([-2 1]);ylim([1.4e4 2.0e4]);
subplot(1,2,2);
plot(squeeze(est_em_n(2,:,:)),...
      lh_em_n(:,:),'b*-'); grid on;
xlabel(pNames{gi(2)});ylabel('lh');
xlim([-7 1]);ylim([1.4e4 2.0e4]);

figure(2); clf;
subplot(1,2,1);
plot(squeeze(est_bfgs_n(1,:,:)),...
      lh_bfgs_n(:,:),'b*-'); grid on;
xlabel(pNames{gi(1)});ylabel('lh');
xlim([-2 1]);ylim([1.4e4 2.0e4]);
subplot(1,2,2);
plot(squeeze(est_bfgs_n(2,:,:)),...
      lh_bfgs_n(:,:),'b*-'); grid on;
xlabel(pNames{gi(2)});ylabel('lh');
xlim([-7 1]);ylim([1.4e4 2.0e4]);



%% simulate

xs = zeros(xDim,N+1);

x = m0;
xs(:,1) = x;
mqw = mean(squeeze(est_bfgs_n(1,end,:)));
mqx = mean(squeeze(est_bfgs_n(2,end,:)));


Q = sinusoid_Q(mqw,mqx);   

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  xs(:,k) = x;
end
endi = find(K<8,1,'last');
figure(1); clf;
subplot(2,1,1);
plot(K(1:endi),H*xs(:,1:endi)); grid on;
subplot(2,1,2);
plot(K(1:endi),Y(:,1:endi)); grid on;

figure(2); clf;

plot(K(1:endi),xs(1,1:endi)/(2*pi)); grid on; title('Freq');



%% Export

%% trajectory

% figure(1); clf;
endi = find(K<10,1,'last');
% plot(K(1:endi),Y(1:endi),'k-x');
% ylim([-0.5 0.5]);
% break;

plt = struct();kw=struct(); kw1 = struct();
kw.alpha = 1.0; kw.ms = 4; kw.mfc = 'black';
kw1.color = '#348ABD';
kw1.lw = 0.9; kw1.alpha = 0.8;
plt.data = {{K(1:endi) Y(1:endi) 'x' kw},{' '}};%,{{K(1:endi) Y(1:endi) 'x' kw}}};
plt.xlabel = '$t$';
plt.ylabel = '$V$';
%plt.legend = {'ECG' ''};
%plt.legendkw = struct('loc','lower center');
plt.w = textwidth*0.7;
plt.margins = [0.0 0.1 0.4 0.35];

%plotstruct(ax,plt);
pyplot('../img/harmonic_trajectory.pdf',plt);




%% Export

textwidth = 426.79134/72.27; % latex textwidth in inches
% Likelihood
plt = struct();kw1=struct();
kw1.color = '#348ABD'; kw1.alpha=0.6; kw1.lw = 0.8;
%(0:itr-1)*avg_time_em;
x = 1:35;
y = lh_em(x,:);

plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5+0.4;%/1.8;
plt.ticklabels = [1 1];
plt.margins = [0.3 0.0 0.4 0.5];

%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
yl = [1e4 2.1e4];
plt.axis = [min(x) max(x) yl];
pyplot('../img/harmonic_em_lh.pdf',plt,'../img/harmonic_em_lh.mat');

plt = struct();
kw1.color = '#E24A33';
x = 1:35;
y = lh_bfgs_n(x,:);
plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5;
plt.ticklabels = [1 0];
plt.margins = [0.3 0.0 0.4 0.1];
%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
plt.axis = [min(x) max(x) yl];
pyplot('../img/harmonic_bf_lh.pdf',plt,'../img/harmonic_bf_lh.mat');


%%

textwidth = 426.79134/72.27; % latex textwidth in inches
% Likelihood
plt = struct();kw1=struct();kw2=struct();
kw1.alpha=0.4; kw1.lw = 0.8;
kw2.alpha=0.4; kw2.ms = 3;
%(0:itr-1)*avg_time_em;
colors = {'#348ABD','#E24A33'};
est = {shiftdim(est_em_n,1),shiftdim(est_bfgs_n,1)};
lh = {lh_em_n,lh_bfgs_n};
xl = [-2 1;
      -7 1];
yl = [1e4 2e4];

plt.ylabel = '$\ell$';
for i=1:2
  y = lh{i};
  y = y+0.02*repmat(0.5-rand(1,size(y,2)),size(y,1),1)*(max(max(y))-min(min(y)));
  xa = est{i};
  
  kw1.color = colors{i};
  kw2.color = colors{i};
  for j=1:2
    x = xa(:,:,j);
    %x = x + 0.01*repmat(0.5-rand(1,size(x,2)),size(x,1),1)*(max(max(x))-min(min(x)));
    plt.data = {{x y '' kw1},{' '}};
    plt.ticklabels = [i-1 2-j];
    plt.margins = [0.2 0.0 i*0.1+(i-1)*0.4 (j-1)*0.1+(2-j)*0.5];
    plt.w = textwidth*0.5+(2-j)*0.4;
    plt.axis = [xl(j,:) yl];
    fn = sprintf('%.0f_%.0f.pdf',i,j);
    fn1 = sprintf('%.0f_%.0f.mat',i,j);
    pyplot(['../img/harmonic_lh_' fn],plt,['../img/harmonic_lh_' fn1]);
  end
end





%% Estimates


plt = struct();
kw1.color = '#348ABD';
kw2.color = '#000000'; kw2.alpha=0.6; kw2.lw = 1.1;
plt.w = textwidth*0.5+0.4;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.50];
x = 1:35;
labelNames = {'\log(\sigma_\omega)', '\log(\sigma_r)','\log(\sigma_x)'};
yl = [min(x) max(x) -1.7 1;
      min(x) max(x) -25 -5;
      min(x) max(x) 0 -7];
for k = 1:numel(gi)
  y1 = squeeze(est_em(k,x,:));
  plt.data = {{x y1 '' kw1},{' '}};
  var = pNames{gi(k)}; 
  plt.ylabel = sprintf('$%s$',labelNames{gi(k)});
  plt.axis = yl(gi(k),:);
  if gi(k) == 3 % lqx is last
    plt.margins(3) = 0.3;
    plt.ticklabels = [1 1];
  else
    plt.margins(3) = 0.1;
    plt.ticklabels = [0 1];
  end
  pyplot(sprintf('../img/harmonic_em_%s.pdf',var),plt,...
         sprintf('../img/harmonic_em_%s.mat',var));
end

plt = struct();
kw1.color = '#E24A33';
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.1];
x = 1:35;
yl = [min(x) max(x) -1.7 1;
      min(x) max(x) -25 -5;
      min(x) max(x) 0 -7];
for k = 1:numel(gi)
  y1 = squeeze(est_bfgs_n(k,x,:));
  plt.data = {{x y1 '' kw1},{' '}};
  var = pNames{gi(k)};
  plt.axis = yl(gi(k),:);
  if gi(k) == 3 % lqx is last
    plt.margins(3) = 0.3;
    plt.ticklabels = [1 0];
  else
    plt.margins(3) = 0.1;
    plt.ticklabels = [0 0];
  end
  %plt.ylabel = sprintf('$%s$',labelNames{k});
  pyplot(sprintf('../img/harmonic_bf_%s.pdf',var),plt,...
         sprintf('../img/harmonic_bf_%s.mat',var));
end

%% Result table
est = {shiftdim(est_em_n,1),shiftdim(est_bfgs_n,1)};
lh = {lh_em_n,lh_bfgs_n};
fem = exp(squeeze(est_em_n(:,end,:)))';
fbf = exp(squeeze(est_bfgs_n(:,end,:)))';
s1 = sqrt(numel(lh{1}(end,:)));
s2 = sqrt(numel(lh{2}(end,:)));

este = [mean(fem)'; 
        std(fem)'/s1; 
        mean(lh{1}(end,:))';
        std(lh{1}(end,:))'/s1];
estb = [mean(fbf)'; 
        std(fbf)'/s2; 
        mean(lh{2}(end,:))';
        std(lh{2}(end,:))'/s2];
filt = lh{2}(end,:)>1.9e4;
s3 = sqrt(numel(lh{2}(end,filt)));

estb2 = [mean(fbf(filt,:))'; 
        std(fbf(filt,:))'/s3; 
        mean(lh{2}(end,filt))';
        std(lh{2}(end,filt))'/s3];


%tr = exp(p_true(gi)); tr(3) = nan;

  rLabels = {'EM($100$)' 'BFGS($100$)' 'BFGS($69$)'};
  cLabels = {'$\sigma_\omega$' '$\pm$' '$\sigma_x$' '$\pm$' '$\ell/\num{e4}$' '$\pm$'};
  
  ord = [1 3 2 4 5 6];
  al = {'S[table-format=1.3,table-omit-exponent]'
        'S[table-format=1.1e2]'
        'S[table-format=1.3,table-omit-exponent]'
        'S[table-format=1.1e2]'
        'S[table-format=1.3,table-omit-exponent,fixed-exponent=4]'
        'S[table-format=1.1e2]'};
  matrix2latex([este(ord)'; estb(ord)'; estb2(ord)'],'heart/Harmonic_results.tex',...
         'alignment',al,'format','%.5e','columnLabels',cLabels,...
         'rowLabels',rLabels,'rowLabelAlignment','r');





