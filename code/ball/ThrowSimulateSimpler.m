%% Setup

global dt P0 H A g0x g0y u

N = 1500;
T = 14;
dt = T/N;





K = (0:N)*dt;
v0 = 300/3.6; % magnitude of the initial velocity
qx = 0.9^2;
qy = qx;%0;%qx/10;

alpha0 = (60/180)*pi; % initial direction
v0x = v0*cos(alpha0);
v0y = v0*sin(alpha0);
g0y = -9.81; % initial y acceleration
g0x = -1; % initial x acceleration

Q = ballisticQ2D(qx,qy);
A1 = [1 dt; 
      0 1];  
A = blkdiag(A1,A1); % for two dimensions
H = zeros(2,4); H(1,1) = 1; H(2,3) = 1;
r = (0.9)^2;
R = r*eye(2);
u = [0 dt*g0x 0 dt*g0y]';

m0 = [0 v0x 0 v0y]';
P0 = eye(size(m0,1));%diag([1e-6 7^-2 1e-6 7^-2]);

%% Simulate

xs = zeros(4,N+1);
%skipY = 4;
ys = zeros(2,N+1);

% simulate
x0 = m0;%mvnrnd(m0,P0)';
atan(x0(4)/x0(2))*180/pi

x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=2:N+1
	x = mvnrnd(A*x,Q)'+u;
	xs(:,k) = x;
	ys(:,k) = mvnrnd(H*x,R)';
end
% the last index where y is still positive
% N = find(xs(4,:)>=0,1,'last')
% K = (0:(N-1))*dt;
% xs = xs(:,1:N);
% ys = ys(:,1:N);
y = ys;
xDim = size(A,1);
MM = zeros(xDim,N+1); MM(:,1) = m0;
PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
dm0 = zeros(xDim,1); dP0 = zeros(xDim);
m = m0; P = P0; lh = 0; 
dm = dm0; dP = dP0;
for k=1:(N+1)
    [m_,P_] = kf_predict(m,P,A,Q,u);

    if k==N+1; break; end; 

    [m,P,~,IM,IS] = kf_update(m_,P_,y(:,k+1),H,R);
    MM(:,k+1) = m;
    PP(:,:,k+1) = P;
    %likelihood(y(:,k+1)-IM,IS)
    lh = lh + likelihood(y(:,k+1)-IM,IS);

end



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





%p0 = {A,Q,H,R};
%[ms,Ps,ms_,Ps_,Ds,lh] = SigmaFilter(p0,ys,[],[],[],[],m0,eye(6));
%[ms,Ss,lh] = SigmaFilterSR(p0,ys,[],[],[],[],m0,eye(6));
 




%% Compute LH and gradient on grid

% parameters are 
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement variance
% set up the starting point
p0 = [v0x v0y qx qy r];
p = p0;
NN = 250;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;


gi = 3;

%as = linspace(0.5*qx,1.5*qx,NN);
as = linspace(0.5*p0(gi),1.5*p0(gi),NN);

for j=1:NN
    %Q = ballisticQ(as(j),qy);
    j
    p(gi) = as(j);
    %p(gi)
    [lh,glh,MM,PP,MM_,PP_] = Ballistic_LH(p,ys,gi);

    lhs(j) = lh;
    glhs(j) = glh;
    
    [MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A); % D = Smoother Gain
    [I1,I2,I3] = EM_I123(A,H,m0,ys,MS,PS,DD);
    glb = EM_LB_Ballistic(p,MS(:,1),gi,N,I1,I2,I3);

    glbs(j) = glb;
end

n = 3;
figure(4); clf;
subplot(n,1,1); 
plot(as,lhs,[p0(gi) p0(gi)],[min(lhs) max(lhs)],'-r'); grid ON; title('Likelihood');
subplot(n,1,2); 
plot(as,glhs,as,glbs,[p0(gi) p0(gi)],[min(glhs) max(glhs)],'-r'); grid ON; title('dLH');
%hold on
%plot(as,glbs,'-g'); grid ON;
subplot(n,1,3); 
plot(as(2:end),diff(lhs)./diff(as)); grid ON; title('numd');

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
%% Test BFGS optimization
NN = 10;
gi = 5;
opt = optimset(@fminunc);
opt.GradObj = 'on';
%opt.TolFun = 1e-15;
opt.TolX = 1e-7;
%opt.Display = 'off';
init = 0.7*qx;
f = @(xx) Ballistic_LH([v0x v0y qx qy xx],ys,gi,-1);
[qqx,VAL,EF,OP,GRAD] = fminunc(f,init,opt);


%% Test EM optimization
gi = 5;
p0 = [v0x v0y qx qy 0.1*r]; % initial guess
p = p0;
%p = {ballisticQ(0.8*qx,qy),R,v0x,v0y};
lh_ = 0;
tol_lh = 1e-6;
tol_delta = 1e-3;
for k=1:1000
  % E-Step
  %lh
  p(gi)
  %pause 
  [lh,~,MM,PP,MM_,PP_] = Ballistic_LH(p,ys);
  if( abs(lh_-lh) < tol_lh && mean(abs(p_(gi)-p(gi))) < tol_delta ); break; end; 
  %lh_-lh
  lh_ = lh;
  p_ = p;
  
  [MS,PS,DD] = rts_smooth2(MM,PP,MM_,PP_,A);
  [I1,I2,I3] = EM_I123(A,H,m0,ys,MS,PS,DD);
  % M-Step
  p = EM_M_Ballistic(p,MS(:,1),gi,N,I1,I2,I3);
  
end
p(gi)



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





