%% r - estimates
funs = plotFuns();
load('../data/Ballistic_r_ux_uy_100_1396.mat');

itr = max_iter_em;
times_em = times_em(1:itr,1:NN);
zerse = times_em <= 0;
% start from zero
times_em = times_em(1:itr,1:NN)-repmat(times_em(1,1:NN),itr,1);
lh_em = lh_em(1:itr,1:NN);
evals_em = cumsum(ones(itr,NN));
evals_em(lh_em <= 0) = 0;
est_em = est_em(:,1:itr,1:NN);
%[lh_em_n,est_em_n] = funs.normalizeBFGS(evals_em,lh_em,est_em);



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
y = lh_em;
plot(y,'b-');
%xlim([0 30]); ylim([18000 24500]);

subplot(2,1,2);
y = lh_bfgs_n;
plot(y,'r');

% EST EM
figure(2); clf;
for k=1:numel(gi)
  subplot(numel(gi),1,k);
  x = 1:max_iter_em;
  y = squeeze(est_em(k,x,:));
  plot(x,y,'-b',x,p_true(gi(k))*ones(size(x)),'-k','lineWidth',2); grid on;
  title(pNames{gi(k)});
end

figure(3); clf;
for k=1:numel(gi)
  subplot(numel(gi),1,k);
  x = 1:max_iter_bfgs;
  y = squeeze(est_bfgs_n(k,x,:));
  plot(x,y,'-r',x,p_true(gi(k))*ones(size(x)),'-k','lineWidth',2); grid on;
  title(pNames{gi(k)});
end
%% trajectory

xs = zeros(4,N+1);
ys = zeros(2,N+1);

% simulate
x = m0;
xs(:,1) = x;
ys(:,1) = H*x;

for k=2:N+1
  x = mvnrnd(A*x,Q)'+u;
  xs(:,k) = x;
	ys(:,k) = mvnrnd(H*x,R)';
  if x(3) < 0; break; end;
end

xs = xs(:,1:k);
ys = ys(:,1:k);
N = k-1;
K = (0:N)*dt;

figure(1); clf;
plot(xs(1,:),xs(3,:),ys(1,:),ys(2,:),'kx');
axis equal; grid on;
%break;

plt = struct();kw1 = struct(); kw = struct();
kw.alpha = 0.5; kw.ms = 4; kw.mfc = 'black'; kw.color = 'black';
kw1.lw = 1.2; kw1.alpha = 0.9; kw1.color = 'white';
plt.data = {{ys(1,:) ys(2,:) 'x' kw},...
            {xs(1,:) xs(3,:) '' kw1}
           };
plt.xlabel = '$\chi$';
plt.ylabel = '$\gamma$';
%plt.legend = {'meas.' 'true'};
%plt.legendkw = struct('loc','lower center');
plt.w = textwidth*0.7;
%plt.h = textwidth*0.7;
plt.margins = [0.1 0 0.35 0.45];

pyplot('../img/ballistic_trajectory.pdf',plt,'../img/ballistic_trajectory.mat');




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
yl = [-2e4 -0.45e4];
plt.axis = [min(x) max(x) yl];
pyplot('../img/ballistic_lh_em.pdf',plt,'../img/ballistic_lh_em.mat');

plt = struct();
kw1.color = '#E24A33';
x = 1:20;%(0:itr-1)*avg_time_bfgs;
y = lh_bfgs_n(x,:);
plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5;
plt.ticklabels = [1 0];
plt.margins = [0.3 0.0 0.4 0.1];
%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
plt.axis = [min(x) max(x) yl];
pyplot('../img/ballistic_lh_bf.pdf',plt,'../img/ballistic_lh_bf.mat');

%%

% Estimates
plt = struct();
kw1.color = '#348ABD';
kw2.color = '#000000'; kw2.alpha=0.6; kw2.lw = 1.1;
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5+0.4;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.50];
x = 1:45;%(0:itr-1)*avg_time_em; 
labelNames = {'', '', '\log(\sigma_r)', '','u_x','u_y'};
yl = [0 0 0 0;0 0 0 0;
      min(x) max(x) -0.6 4;
      0 0 0 0;
      min(x) max(x) 0.2 -3.5;
      min(x) max(x) 0.2 -20];
for k = 1:3
  y1 = squeeze(est_em(k,x,:));
  plt.data = {{x y1 '' kw1},{x p_true(gi(k))*ones(size(x)) '' kw2}};
  var = pNames{gi(k)};
  plt.ylabel = sprintf('$%s$',labelNames{gi(k)});
  plt.axis = yl(gi(k),:);
  if gi(k) == 3 % r is last
    plt.margins(3) = 0.2;
    plt.ticklabels = [1 1];
  else
    plt.margins(3) = 0.1;
    plt.ticklabels = [0 1];
  end
  pyplot(sprintf('../img/ballistic_em_%s.pdf',var),plt,...
         sprintf('../img/ballistic_em_%s.mat',var));
end

plt = struct();
kw1.color = '#E24A33';
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.1];
x = 1:30;%1:max_iter_bfgs;%(0:itr-1)*avg_time_em; 
yl = [0 0 0 0;0 0 0 0;
      min(x) max(x) -0.6 4;
      0 0 0 0;
      min(x) max(x) 0.2 -3.5;
      min(x) max(x) 0.2 -20];
for k = 1:3
  y1 = squeeze(est_bfgs_n(k,x,:));
  plt.data = {{x y1 '' kw1},{x p_true(gi(k))*ones(size(x)) '' kw2}};
  var = pNames{gi(k)};
  plt.axis = yl(gi(k),:);
  if gi(k) == 3 % r is last
    plt.margins(3) = 0.2;
    plt.ticklabels = [1 0];
  else
    plt.margins(3) = 0.1;
    plt.ticklabels = [0 0];
  end
  %plt.ylabel = sprintf('$%s$',labelNames{k});
  pyplot(sprintf('../img/ballistic_bf_%s.pdf',var),plt,...
         sprintf('../img/ballistic_bf_%s.mat',var));
end

%% Result table
este = mean(squeeze(est_em(:,end,:)),2); este(1) = exp(este(1));
este(4) = mean(lh_em(end,:));
estb = mean(squeeze(est_bfgs_n(:,end,:)),2); estb(1) = exp(estb(1));
estb(4) = mean(lh_bfgs_n(end,:));
tr = p_true(gi); tr(1) = exp(tr(1)); tr(4) = nan;

  rLabels = {'BFGS' 'EM' 'True'};
  cLabels = {'$g_\chi$' '$g_\gamma$' '$\sigma_r$' '$\ell/\num{e3}$'};
  
  ord = [2 3 1 4];
  al = {'S' 'S' 'S' 'S[fixed-exponent=3]'};
  matrix2latex([estb(ord)'; este(ord)'; tr(ord);],'ball/Ballistic_results.tex',...
         'alignment',al,'format','%.5e','columnLabels',cLabels,...
         'rowLabels',rLabels,'rowLabelAlignment','r');
       
%%

rn = rand(1,100);
p0 = p_true;
plot(1:100,p0(3)*(2-2*rn),1:100,p0(3)*ones(size(rn)));














