% test both ways of computing 
% the likelihood and its gradient

%% AR(1) model

N = 150;
A = 1;
Q = 1^2;
R = Q*2.2;
H = 1;

K = 0:N;
m0 = 0;
P0 = 1;
%% Simulate

xs = zeros(1,N+1);
ys = zeros(1,N+1);

% simulate
x0 = m0;
x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=2:N+1
	x = A*x+sqrt(Q)*randn;
  xs(:,k) = x;
	ys(:,k) = x+sqrt(R)*randn;
end


figure(1); clf;
plot(K,xs,K,ys,'kx');

%% compute likelihood for A on grid


NN = 250;
alpha = 2.0;
as = linspace((1-alpha)*A,(1+alpha)*A,NN); % for q1
lhs = zeros(1,NN);
y = ys;
m = m0; P = P0;
for j=1:NN
  lh = 0;
  for k=1:(N+1)

      [m_,P_] = kf_predict(m,P,as(j),Q);

      if k==N+1; break; end; 

      [m,P,~,my,Sy] = kf_update(m_,P_,y(:,k+1),H,R);

      lh = lh + likelihood(y(:,k+1)-my,Sy);

  end
  lhs(j) = lh;
end

figure(2); clf;
plot(as,lhs); grid on;
yl = ylim;

%% demonstrate EM

KK = 2;
jj = floor(NN/10);
A_ = as(jj);
lbs = zeros(KK,NN);
jjs = zeros(1,KK+1);
jjs(1) = jj;
xDim = size(A,1);
for kk=1:KK
  MM = zeros(xDim,N+1); MM(:,1) = m0;
  PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
  MM_ = zeros(xDim,N+1);
  PP_ = zeros(xDim,xDim,N+1);
  m = m0; P = P0;
  for k=1:(N+1)
      [m_,P_] = kf_predict(m,P,A_,Q);
      MM_(:,k) = m_;
      PP_(:,:,k) = P_;

      if k==N+1; break; end; 

      [m,P,~,my,Sy] = kf_update(m_,P_,y(:,k+1),H,R);
      MM(:,k+1) = m;
      PP(:,:,k+1) = P;
  end
  [MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A_); % D = Smoother Gain
  
  % evaluate lower bound
  for j=1:NN
    [I1,I2,I3] = EM_I123(as(j),H,m0,ys,MS,PS,DD);
    I1 =   log(det(P0))+trace(I1/P0);
    I2 = N*log(det(Q))+trace(I2/Q);
    I3 = N*log(det(R))+trace(I3/R);
    lbs(kk,j) = -0.5*(I1+I2+I3);
    if j==jj
      df = lhs(j)-lbs(kk,j);
    end
  end
  lbs(kk,:) = lbs(kk,:)+df;
  [~,jj] = max(lbs(kk,:));
  A_ = as(jj);
  jjs(kk+1) = jj;
end
figure(3); clf;
plot(as,lhs,as,lbs,'--r',as(jjs),lhs(jjs),'k*'); ylim(yl); grid on;


%% export plots

textwidth = 426.79134/72.27; % latex textwidth in inches
% plot the true locations and the measurements


plt = struct();kw=struct();
%kw.alpha = 0.2; kw.linewidth = 0.5;
kw.ms = 3; kw.alpha = 0.8; 
plt.data = {{K xs},{K ys 'kx' kw}};
plt.w = textwidth*0.5;%/1.8;
plt.xlabel = '$k$';
plt.legend = {'$x$' '$y$'};
%plt.title = '(a) AR(1) simulation';
plt.alpha = 0.1;
pyplot('../img/ar1_ex_a.pdf',plt,'../img/ar1_ex_a.mat');

plt = struct();kw=struct();
plt.data = {{as lhs},{''}};
plt.w = textwidth*0.5;
plt.xlabel = '$\theta$';
%plt.legend = {'$x$' '$y$'};
%plt.title = '(b) Likelihood of $a$';
plt.alpha = diag([0.05 0.05 0.15 0.15]);
pyplot('../img/ar1_ex_b.pdf',plt,'../img/ar1_ex_b.mat');

lbsmax = zeros(1,KK);
lbsmax(1) = lhs(jjs(1));
for k=1:KK
  lbsmax(k+1) = lbs(k,jjs(k+1));
end
plt = struct();kw=struct();
kw2.color = '#E24A33';
plt.data = {{as lhs},{as lbs' '--' kw2},{as(jjs),lbsmax,'k*'}};
plt.w = textwidth*0.7;
plt.xlabel = '$\theta$';
plt.legend = {'$\ell$' '$\mathcal{B}$'};
plt.alpha = diag([0.05 0.05 0.15 0.15]);
pyplot('../img/ar1_ex_em.pdf',plt,'../img/ar1_ex_em.mat');








