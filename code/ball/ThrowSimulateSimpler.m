%% Setup

global dt m0 P0 H A


N = 500;
T = 14;
dt = T/N;


%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%
lqx = 0.8;       % parameterization in logarithm of standard deviation
lqy = lqx/2;
lr =  log(2.5);


K = (0:N)*dt;
v0 = 300/3.6; % magnitude of the initial velocity


alpha0 = (60/180)*pi; % initial direction
v0x = v0*cos(alpha0);
v0y = v0*sin(alpha0);
g0y = -9.81; % initial y acceleration
g0x = -1.5; % initial x acceleration
A1 = [1 dt; 
      0 1];  
A = blkdiag(A1,A1); % for two dimensions
H = zeros(2,4); H(1,1) = 1; H(2,3) = 1;



Q = ballisticQ2D(lqx,lqy);
SQ = chol(Q,'lower');
h = @(x) H*x;
R = ballisticR(lr);
SR = chol(R,'lower');
u = ballisticU(g0x,g0y);


m0 = [0 v0x 0 v0y]';
P0 = eye(size(m0,1));%diag([1e-6 7^-2 1e-6 7^-2]);

%u = zeros(size(m0,1),1);
f = @(x) A*x+repmat(u,1,size(x,2));


%% Simulate
[usig,w] = CKFPoints(size(A,1));
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
xs = zeros(4,N+1);
%skipY = 4;
ys = zeros(2,N+1);

% simulate
x0 = m0;%mvnrnd(m0,P0)';
%atan(x0(4)/x0(2))*180/pi

x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=2:N+1
	x = mvnrnd(A*x,Q)'+u;
  xs(:,k) = x;
	ys(:,k) = mvnrnd(H*x,R)';
end
%u = zeros(size(x));
y = ys;
xDim = size(A,1);
MM = zeros(xDim,N+1); MM(:,1) = m0;
PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
MM_ = zeros(xDim,N+1);
PP_ = zeros(xDim,xDim,N+1);

dm0 = zeros(xDim,1); dP0 = zeros(xDim);
m = m0; P = P0; lh = 0; m2 = m0; S = P0;
dm = dm0; dP = dP0;
for k=1:(N+1)
    [m_,P_] = kf_predict(m,P,A,Q,[],u);
    [m2_,S_] = SigmaKF_Predict(m2,S,f,SQ,usig,w);
    MM_(:,k) = m_;
    PP_(:,:,k) = P_;
    
    if k==N+1; break; end; 

    [m,P,~,my,Sy] = kf_update(m_,P_,y(:,k+1),H,R);
    [m2,S] = SigmaKF_Update(m2_,S_,y(:,k+1),h,SR,usig,w);
    
    
    if k < 10
      P2_ = S_*S_';
      P2 = S*S';
      disp([rmse(m_,m2_)
            rmse(P_(:),P2_(:))
            rmse(m,m2)
            rmse(P(:),P2(:))]);
      
    end
    
    MM(:,k+1) = m;
    PP(:,:,k+1) = P;
    %likelihood(y(:,k+1)-IM,IS)
    %lh = lh + likelihood(y(:,k+1)-IM,IS);

end
[MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A); % D = Smoother Gain



figure(1); clf;
plot(xs(1,:),xs(3,:),ys(1,:),ys(2,:),'kx',MM(1,:),MM(3,:));
axis equal; grid on;
figure(2); clf;
subplot(2,2,1);
plot(K,xs(1,:)); 
subplot(2,2,2);
plot(K,xs(3,:)); 
subplot(2,2,3);
plot(K,xs(2,:));
subplot(2,2,4);
plot(K,xs(4,:));

figure(3); clf; n=3;
subplot(n,1,1);
plot(K,squeeze(PP(1,1,:)));
subplot(n,1,2);
plot(K,squeeze(PP(2,2,:)));
subplot(n,1,3);
plot(K,squeeze(PP(3,3,:)));

figure(4); clf;
plot(K,sqrt(mean((xs-MM).^2)),K,sqrt(mean((xs-MS).^2)));



%p0 = {A,Q,H,R};
%[ms,Ps,ms_,Ps_,Ds,lh] = SigmaFilter(p0,ys,[],[],[],[],m0,eye(6));
%[ms,Ss,lh] = SigmaFilterSR(p0,ys,[],[],[],[],m0,eye(6));
 




%% Compute LH and gradient on grid

% parameters are 
%
% p{1}=lqx,  x process variance, or shared if gi = 4
% p{2}=lqy,  y process variance
% p{3}=lr,   measurement variance
% p{4}=lqxqy,   joint process variance
% p{5}=ux,    constant input x-component
% p{6}=uy,    constant input y-component 
p0 = [lqx lqy lr 0 g0x g0y];
p = p0;
NN = 20;
gi = [3 5];
alpha = 0.2;
as = zeros(numel(gi),NN^(numel(gi)));
%as = linspace(0.5*qx,1.5*qx,NN);
for k = 1:numel(gi)

  true = p0(gi(k));
  if gi(k) < 5
    as1 = linspace(true+log(0.5),true-log(0.5),NN);
  else
    as1 = linspace((1-alpha)*true,(1+alpha)*true,NN);
  end
  if numel(gi) > 1
    if k == 1
      as(k,:) = kron(ones(size(as1)),as1);
    else
      as(k,:) = kron(as1,ones(size(as1)));
    end
  else
    as = as1;
  end
end

lhs = zeros(1,size(as,2)); 
glhs = zeros(size(as)); 
glbs = glhs;



for j=1:size(as,2)
    j
    for k=1:numel(gi)
        p(gi(k)) = as(k,j);
    end

    [lh,glh,MM,PP,MM_,PP_,Q,u] = Ballistic_LH(p,ys,gi);
    lhs(1,j) = lh;
    glhs(:,j) = glh;
    %glh(2)
    
    [MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A); % D = Smoother Gain
    [I1,I2,I3] = EM_I123(A,H,m0,ys,MS,PS,DD,u);
    smk = sum(MS(:,2:end),2);
    smkk = sum(MS(:,1:end-1),2);
    glb = EM_LB_Ballistic(p,MS(:,1),gi,N,I1,I2,I3,smk,smkk,Q,u);
    glb(1)-glh(1)
    glb(2)-glh(2)
    glbs(:,j) = glb;
    % Sigma Filter and Smoother
%     [lh,glh,MMM,SSS] = Ballistic_LH_Sigma(p,ys,gi);
%     lhs(1,j) = lh;
%     glhs(1,j) = glh;
%     SQ = chol(Q,'lower');
%     [MMS,SS,DD] = SigmaSmoothSR(MMM,SSS,f,SQ,usig,w); % D = Smoother Gain
%     [I1,I2,I3] = EM_I123_Sigma(f,h,m0,ys,MMS,SS,DD);
%     smk = sum(MMS(:,2:end),2);
%     smkk = sum(MMS(:,1:end-1),2);
%     glbs(1,j) = EM_LB_Ballistic(p,MS(:,1),gi,N,I1,I2,I3,smk,smkk,Q);

    % SR filter and RTS smoother
%     [lh,glh,MM,SS,MM_,SS_] = Ballistic_LH_SR(p,ys,gi);
%     lhs(1,j) = lh;
%     glhs(1,j) = glh;
%     [MS,PS,DD] = rts_smooth_sr(MM,SS,MM_,SS_,A); % D = Smoother Gain
%     [I1,I2,I3] = EM_I123(A,H,m0,ys,MS,PS,DD,u);
%     glbs(1,j) = EM_LB_Ballistic(p,MS(:,1),gi,N,I1,I2,I3);
    % Sigma Filter and Smoother

    
end

if numel(gi) == 1
  n = 2; m= 1;
  figure(1); clf;
  subplot(n,m,1);
  if gi < 5
    eas = exp(as);
  else
    eas = as;
  end
  plot(eas,lhs'); grid on;
  subplot(n,m,2);
  diff1 = diff(lhs(1,:))./diff(as);
  plot(eas,glhs,eas,glbs,eas(2:end),diff1); grid on;
  legend('sens','em','num');
else
  X = reshape(as(1,:),NN,NN);
  Y = reshape(as(2,:),NN,NN);
  Z = reshape(lhs,NN,NN);
  dlh1 = reshape(glhs(1,:),NN,NN);
  dlh2 = reshape(glhs(2,:),NN,NN);
  dlh = dlh1.*X+dlh2.*Y;
  
  dem1 = reshape(glbs(1,:),NN,NN);
  dem2 = reshape(glbs(2,:),NN,NN);
  dem = dem1.*X+dem2.*Y;
  
  
  figure(1); clf;
  plot(exp(X(:,1)),dem1(:,1),exp(X(:,1)),dlh1(:,1))
  figure(2); clf;
  plot(Y(1,:),dem1(1,:),Y(1,:),dlh1(1,:))
  figure(3); clf;
  plot(exp(X(:,1)),dem2(:,1),exp(X(:,1)),dlh2(:,1))
  figure(4); clf;
  plot(Y(1,:),dem2(1,:),Y(1,:),dlh2(1,:))
  
end
%figure(2); clf;
% plot(eas,glhs-glbs);

% n = 4; m= 2;
% figure(1); clf;
% subplot(n,m,1);
% eas = exp(as);
% plot(eas,lhs'); grid on;
% subplot(n,m,3);
% plot(eas,glhs'); grid on;
% subplot(n,m,5);
% plot(eas,glbs'); grid on;
% subplot(n,m,2);
% %plot(K,MM(1,:)-MMM(1,:));
% plot(eas,sqrt((lhs(1,:)-lhs(2,:)).^2));
% subplot(n,m,4);
% %plot(K,squeeze(PP(1,1,:))-squeeze(SSS(1,1,:)).^2);
% plot(eas,sqrt((glhs(1,:)-glhs(2,:)).^2));
% subplot(n,m,6);
% %plot(K,squeeze(PS(1,1,:))-squeeze(SS(1,1,:)).^2,K,MS(1,:)-MMS(1,:));
% plot(eas,sqrt((glbs(1,:)-glbs(2,:)).^2));
% subplot(n,m,7);
%plot(K,squeeze(PS(1,1,:))-squeeze(SS(1,1,:)).^2,K,MS(1,:)-MMS(1,:));
%diff1 = diff(lhs(1,:))./diff(as);
% diff2 = diff(lhs(2,:))./diff(as);
% plot(eas(1:end-1),diff1,eas(1:end-1),diff2); grid on;
% subplot(n,m,8);
% plot(eas(1:end-1),sqrt((diff2-glhs(2,2:end)).^2),eas(1:end-1),sqrt((diff2-glbs(2,2:end)).^2),...
%   eas(1:end-1),sqrt((diff1-glhs(1,2:end)).^2),eas(1:end-1),sqrt((diff1-glbs(1,2:end)).^2)); grid on;


% figure(2); clf;
% subplot(2,1,1);
% plot(eas(1:end-1),sqrt((diff1-glhs(1,2:end)).^2),eas(1:end-1),sqrt((diff1-glbs(1,2:end)).^2)); grid on;
% title('KF/RTS');
% subplot(2,1,2);
% plot(eas(1:end-1),sqrt((diff2-glhs(2,2:end)).^2),eas(1:end-1),sqrt((diff2-glbs(2,2:end)).^2)); grid on;
% title('Sigma');




%plot(glbs(2,:));

%plot(as,lhs,[p0(gi) p0(gi)],[min(lhs) max(lhs)],'-r'); grid ON; title('Likelihood');
%subplot(n,1,2); 
%plot(as,glhs,as,glbs,[p0(gi) p0(gi)],[min(glhs) max(glhs)],'-r'); grid ON; title('dLH');
%hold on
%plot(as,glbs,'-g'); grid ON;
%subplot(n,1,3); 
%plot(as(2:end),diff(lhs)./diff(as)); grid ON; title('numd');

%figure(2);clf;

% figure(3); clf;
% subplot(3,1,1); 
% plot(as,lhsSR); grid ON; title('LikelihoodSR');
% subplot(3,1,2); 
% plot(as,lbsSR); grid ON; title('Lower bound SR');
% subplot(3,1,3);
% plot(as,glbsSR); grid ON; title('d Lower bound SR');
% figure(4);clf;
% plot(as,lhsSR-lbsSR);grid on;
%% Test optimizations
gi = 5;
p_true = [v0x v0y lqx lqy r]; % initial guess
true = p_true(gi);
min_iter_em =   10;
max_iter_em =   10;
min_iter_bfgs = 15;
max_iter_bfgs = 15;
NN = 20;
est_em =   zeros(numel(gi),max_iter_em,NN);
est_bfgs = zeros(numel(gi),max_iter_bfgs,NN);

evals_em = zeros(1,NN);
evals_bfgs = zeros(2,NN);

xs = zeros(4,N+1);
ys = zeros(2,N+1);

% simulate
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
  
  % INITIAL POINT
  p0 = p_true;
  p0(gi) = p0(gi)*(rand+0.5);
  
  % EM
  tic;
  [~,~,vals] = Ballistic_EM(p0,gi,ys,[],[],max_iter_em,min_iter_em);
  tm = toc;
  est_em(:,:,k) = vals;
  evals_em(1,k) = tm;
  
  % BFGS
  tic;
  [~,~,vals,fcn_evals] = Ballistic_BFGS(p0,gi,ys,[],[],max_iter_bfgs,min_iter_bfgs);
  tm = toc;
  num = size(vals,2);
  est_bfgs(:,1:num,k) = vals;
  if num < max_iter_bfgs
    est_bfgs(:,num+1:end,k) = repmat(vals(:,end),1,max_iter_bfgs-num);
  end
  evals_bfgs(:,k) = [tm;fcn_evals];
end



figure(1); clf;
est_em1 =   reshape(est_em,[max_iter_em NN  1]);
est_bfgs1 = reshape(est_bfgs,[max_iter_bfgs NN 1]);
subplot(2,1,1);
plot(1:max_iter_em,-1*(true-est_em1)./true,'-b');
ylim([-0.5,0.5]);xlim([1,max_iter_em]);
subplot(2,1,2);
plot(1:max_iter_bfgs,(-1*(true-est_bfgs1)./true),'-b');
ylim([-0.5,0.5]);xlim([1,max_iter_bfgs]);

%% Plot

textwidth = 426.79134/72.27; % latex textwidth in inches
% plot the true locations and the measurements
plt = struct();kw=struct();
kw.alpha = 0.9;
plt.data = {{xs(1,:) xs(4,:) '' kw},...
			{ms(1,:) ms(4,:) '' kw},...
			{mF(1,:) mF(4,:) '' kw},...
			{ys(1,yI) ys(2,yI) 'x'},...
			};
plt.xlabel = '$\mathrm{[m]}$';
plt.ylabel = '$\mathrm{[m]}$';
plt.legend = {'true' 'filter' 'smoother' 'measurement'};
plt.legendkw = struct('loc','lower center');
plt.w = textwidth*0.5;
%plotstruct(ax,plt);
pyplot('../img/ex1_pos_meas.pdf',plt);

plt = struct();
fmean = xs(2,:)-ms(2,:);
smean = xs(2,:)-mF(2,:);
ferr = 2*sqrt(squeeze(Ps(2,2,:)));
serr = 2*sqrt(squeeze(PF(2,2,:)));
plt.ylabel = '$\dot{x}_{\mathrm{true}}-\dot{x}_{\mathrm{mean}}\,\mathrm{[m]}$';
plt.xlabel = '$t\,\mathrm{[s]}$';
plt.legend = {'$\mathrm{Err}_f$' '$\mathrm{Err}_s$'};
plt.legendkw = struct('loc','upper right');
plt.data = {{K fmean '' struct('yerr',ferr)},{K smean '' struct('yerr',serr)}};
plt.w = textwidth*0.5;
plotstruct(plt);
pyplot('../img/ex1_err.pdf',plt)

%plt.data = {{K xs(5,:)},{}};
%plt.ylabel = '$\dot{y}$';
%pyplot('../img/ex1_x5.pdf',plt);
%plt.data = {{K xs(6,:)},{}};
%plt.ylabel = '$\ddot{y}$';
%pyplot('../img/ex1_x6.pdf',plt);

break

plt = struct()
plt.x = [ms(1,:)' mF(1,:)'];
plt.y = [ms(4,:)' mF(4,:)'];
plt.xlabel = '$x$';
plt.ylabel = '$y$';
plt.legend = {'filtered' 'smoothed'};
plt.w = 4;

pyplot('../img/test.pdf',plt);
error('stop');

figure;
plot(xs(1,:),xs(4,:),ys(1,yI),ys(2,yI),'xk');
figure;
boundedline(ms(1,:),ms(4,:),[1.97*sqrt(squeeze(Ps(1,1,:))) 1.97*sqrt(squeeze(Ps(4,4,:)))])
figure;
boundedline(mF(1,:),mF(4,:),[1.97*sqrt(squeeze(PF(1,1,:))) 1.97*sqrt(squeeze(PF(4,4,:)))])


figure;
subplot(3,2,1);
plot(K,xs(1,:)); hold on;
boundedline(K,ms(1,:),1.97*sqrt(squeeze(Ps(1,1,:))));
plot(K(yI),ys(1,yI),'kx');
subplot(3,2,2);
plot(K,xs(4,:)); hold on;
boundedline(K,ms(4,:),1.97*sqrt(squeeze(Ps(4,4,:))));
plot(K(yI),ys(2,yI),'kx');
subplot(3,2,3);
plot(K,xs(2,:));
subplot(3,2,4);
plot(K,xs(5,:));
subplot(3,2,5);
plot(K,xs(3,:));
subplot(3,2,6);
plot(K,xs(6,:));

figure;
subplot(2,2,1);
plot(K,squeeze(Ps(1,1,:)));
subplot(2,2,2);
plot(K,squeeze(Ps(4,4,:)));
subplot(2,2,3);
plot(K,squeeze(Ps(2,2,:)));
subplot(2,2,4);
plot(K,squeeze(Ps(5,5,:)));

figure;
subplot(2,2,1);
plot(K,squeeze(PF(1,1,:)));
subplot(2,2,2);
plot(K,squeeze(PF(4,4,:)));
subplot(2,2,3);
plot(K,squeeze(PF(2,2,:)));
subplot(2,2,4);
plot(K,squeeze(PF(5,5,:)));





