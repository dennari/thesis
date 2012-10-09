%% r - estimates
funs = plotFuns();
load('../data/Ballistic_r_ux_uy_20_1500.mat');

itr = max_iter_em-1;
times_em = times_em(1:itr,:);
% start from zero
times_em = times_em-repmat(times_em(1,:),itr,1);
lh_em = lh_em(1:itr,:);
est_em = est_em(:,1:itr,:);

times_bfgs = times_bfgs(2:itr+1,:);
zers = times_bfgs <= 0;
times_bfgs = times_bfgs-repmat(times_bfgs(1,:),itr,1);
lh_bfgs = lh_bfgs(2:itr+1,:);
est_bfgs = est_bfgs(:,2:itr+1,:);
evals_bfgs = evals_bfgs(2:itr+1,:);


avg_time_em = mean(mean(diff(times_em(1:itr,:),1,1)));
df=diff(times_bfgs,1,1);
avg_time_bfgs = mean(df(df>1.1&df<1.2)); % choose a good interval

[lh_bfgs_n,est_bfgs_n] = funs.normalizeBFGS(evals_bfgs,lh_bfgs,est_bfgs);
%break
% LH
figure(1); clf;
% normalized
subplot(3,1,1);
plot((0:itr-1)*avg_time_em,lh_em,'-b',(0:itr-1)*avg_time_em,lh_bfgs_n,'-r');
xlim([0 15]); ylim([-0.5e5 0.2e5]);

% original
subplot(3,1,2);
times_bfgs(zers) = nan; lh_bfgs(zers) = nan;
plot(times_em,lh_em,'-b',times_bfgs,lh_bfgs,'-r');
xlim([0 35]); ylim([-1e5 8000]);

% mean over runs
subplot(3,1,3);
plot((0:itr-1)*avg_time_em,mean(lh_em,2),'-b',(0:itr-1)*avg_time_bfgs,mean(lh_bfgs_n,2),'-r');
xlim([0 35]); ylim([-1e5 8000]);

% EST EM
figure(2); clf;
subplot(3,1,1);
x = (0:itr-1)*avg_time_em; 
y1 = squeeze(est_em(1,:,:));
plot(x,y1,'-b'); grid on;
subplot(3,1,2);
y1 = squeeze(est_em(2,:,:));
plot(x,y1,'-b'); grid on;
subplot(3,1,3);
y1 = squeeze(est_em(3,:,:));
plot(x,y1,'-b'); grid on;


% EST BF
figure(3); clf;
subplot(3,1,1);
x = (0:itr-1)*avg_time_bfgs; 
y1 = squeeze(est_bfgs_n(1,:,:));
plot(x,y1,'-b'); grid on;
subplot(3,1,2);
y1 = squeeze(est_bfgs_n(2,:,:));
plot(x,y1,'-b'); grid on;
subplot(3,1,3);
y1 = squeeze(est_bfgs_n(3,:,:));
plot(x,y1,'-b'); grid on;

%% trajectory
global dt
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
R = ballisticR(r);
u = ballisticU(g0x,g0y);
m0 = [0 v0x 0 v0y]';
P0 = eye(size(m0,1));

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
kw1.lw = 1; kw1.alpha = 0.7; kw1.color = 'white';
plt.data = {{ys(1,:) ys(2,:) 'x' kw},...
            {xs(1,:) xs(3,:) '' kw1}
           };
plt.xlabel = '$\chi$';
plt.ylabel = '$\gamma$';
%plt.legend = {'meas.' 'true'};
%plt.legendkw = struct('loc','lower center');
plt.w = textwidth*0.7;
%plt.h = textwidth*0.7;
plt.margins = [0 0 0.35 0.35];

pyplot('../img/ballistic_trajectory.pdf',plt,'../img/ballistic_trajectory.mat');




%% Export

textwidth = 426.79134/72.27; % latex textwidth in inches
% plot the true locations and the measurements
%%
% Likelihood
plt = struct();kw1=struct();
kw1.color = '#348ABD'; kw1.alpha=0.8; kw1.lw = 1.2;
x = (0:itr-1)*avg_time_em;
y = lh_em;
plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5+0.4;%/1.8;
plt.ticklabels = [0 1];
plt.margins = [0.0 0.0 0.1 0.5];
%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
yl = [-0.5e5 0.2e5];
plt.axis = [min(x) max(x) yl];
pyplot('../img/ballistic_lh_em.pdf',plt,'../img/ballistic_lh_em.mat');

plt = struct();kw1=struct();
kw1.color = '#E24A33'; kw1.alpha=0.8; kw1.lw = 1.2;
x = (0:itr-1)*avg_time_bfgs;
y = lh_bfgs_n;
plt.data = {{x y '' kw1},{' '}};
plt.w = textwidth*0.5;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.1];
%plt.xlabel = '$k$'; 
plt.ylabel = '$\ell$';
%plt.alpha = 0.1;
plt.axis = [min(x) max(x) yl];
pyplot('../img/ballistic_lh_bf.pdf',plt,'../img/ballistic_lh_bf.mat');


%% Estimates
plt = struct();kw1=struct();
kw1.color = '#348ABD'; kw1.alpha=0.8; kw1.lw = 1.2;
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5+0.4;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.50];
x = (0:itr-1)*avg_time_em; 
labelNames = {'', '', '\log(\sigma_r)', '','u_x','u_y'};
yl = [0 0 0 0;0 0 0 0;
      min(x) max(x) -5 8;
      0 0 0 0;
      min(x) max(x) -2 2;
      min(x) max(x) -20 20];
for k = 1:3
  y1 = squeeze(est_em(k,:,:));
  plt.data = {{x y1 '' kw1},{' '}};
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

plt = struct();kw1=struct();
kw1.color = '#E24A33'; kw1.alpha=0.8; kw1.lw = 1.2;
%plt.xlabel = '$k$'; 
plt.w = textwidth*0.5;
plt.ticklabels = [0 0];
plt.margins = [0.0 0.0 0.1 0.1];
x = (0:itr-1)*avg_time_em; 
yl = [0 0 0 0;0 0 0 0;
      min(x) max(x) -5 8;
      0 0 0 0;
      min(x) max(x) -2 2;
      min(x) max(x) -20 20];
for k = 1:3
  y1 = squeeze(est_bfgs_n(k,:,:));
  plt.data = {{x y1 '' kw1},{' '}};
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





