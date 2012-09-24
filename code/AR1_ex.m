% test both ways of computing 
% the likelihood and its gradient

%% AR(1) model

N = 250;
A = 1;
Q = 1^2;
R = Q/2^2;
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



% compute likelihood for A on grid


NN = 150;
alpha = 0.5;
as = linspace((1-alpha)*A,(1+alpha)*A,NN); % for q1
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;
y = ys;
xDim = size(A,1);
MM = zeros(xDim,N+1); MM(:,1) = m0;
PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
MM_ = zeros(xDim,N+1);
PP_ = zeros(xDim,xDim,N+1);
m = m0; P = P0;
for j=1:NN
  lh = 0;
  for k=1:(N+1)

      [m_,P_] = kf_predict(m,P,as(j),Q);
      MM_(:,k) = m_;
      PP_(:,:,k) = P_;

      if k==N+1; break; end; 

      [m,P,~,my,Sy] = kf_update(m_,P_,y(:,k+1),H,R);


      MM(:,k+1) = m;
      PP(:,:,k+1) = P;
      lh = lh + likelihood(y(:,k+1)-my,Sy);

  end
  lhs(j) = lh;
end

figure(2); clf;
plot(as,lhs); grid on;


%% export plots

textwidth = 426.79134/72.27; % latex textwidth in inches
% plot the true locations and the measurements


plt = struct();kw=struct();
%kw.alpha = 0.2; kw.linewidth = 0.5;
plt.data = {{K xs},{K ys 'kx'}};
plt.w = textwidth*0.5;
plt.A = A;
pyplot('../img/ar1_ex.mat',plt);










