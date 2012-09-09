%% Setup

global dt P0 H A g0x g0y

T = 15;
N = 1000;
dt = T/N;
K = (0:N)*dt;

qx = (7*dt)^2;
qy = (dt)^2;
v0 = 40; % magnitude of the initial velocity
alpha0 = 50/360*2*pi; % initial direction
v0x = v0*cos(alpha0);
v0y = v0*sin(alpha0);
g0y = -9.81; % initial y acceleration
g0x = 0; % initial x acceleration

Q = ballisticQ(qx,qy);
A1 = [1 dt dt^2/2; 0 1 dt; 0 0 1]; % for one dimension
A = blkdiag(A1,A1); % for two dimensions
H = zeros(2,6); H(1,1) = 1; H(2,4) = 1;
r = 2.5;
R = r^2*eye(2);


m0 = [0 v0x g0x 0 v0y g0y]';
P0 = diag([1e-6 7^2 1e-6 1e-6 7^2 1e-6]);

%% Simulate

xs = zeros(6,N+1);
%skipY = 4;
ys = zeros(2,N+1);

% simulate
x0 = mvnrnd(m0,P0)';
norm([x0(5) x0(2)])
atan(x0(5)/x0(2))*360/(2*pi)

x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=1:N
	x = mvnrnd(A*x,Q)';
	xs(:,k+1) = x;
	%ys(:,k+1) = NaN;
	%if ~mod(k,skipY)
	ys(:,k+1) = mvnrnd(H*x,R)';
	%end
end
% the last index where y is still positive
N = find(xs(4,:)>=0,1,'last');
K = (0:(N-1))*dt;
xs = xs(:,1:N);
ys = ys(:,1:N);
%yI = ~isnan(ys(1,:));

figure(1); clf;
plot(xs(1,:),xs(4,:),ys(1,:),ys(2,:),'kx');
figure(2); clf;
subplot(3,2,1);
plot(K,xs(1,:)); 
subplot(3,2,2);
plot(K,xs(4,:)); 
subplot(3,2,3);
plot(K,xs(2,:));
subplot(3,2,4);
plot(K,xs(5,:));
subplot(3,2,5);
plot(K,xs(3,:));
subplot(3,2,6);
plot(K,xs(6,:));
%p0 = {A,Q,H,R};
%[ms,Ps,ms_,Ps_,Ds,lh] = SigmaFilter(p0,ys,[],[],[],[],m0,eye(6));
%[ms,Ss,lh] = SigmaFilterSR(p0,ys,[],[],[],[],m0,eye(6));


%% Compute LH and gradient on grid

% parameters are 
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx, x process variance
% p{4}=qy, y process variance
% p{5}=r, measurement variance
% set up the starting point
p{1} = v0x;
p{2} = v0y;
p{3} = qx;
p{4} = qy;
p{5} = r;

h = @(x,k,p) H*x;
f = @(x,k,p) A*x;

NN = 80;
as = linspace(2,3,NN);
lhs = zeros(1,NN);
lhsSR = lhs;
glhs = lhs;
lbs = lhs;
lbsSR = lhs;
glbs = lhs;
glbsSR = lhs;
p0 = {A,Q,H,R};
for k=1:NN
    %k
    %as(k)
    %p0{2} = ballisticQ(as(k),qy);
    p0{4} = as(k)*eye(2);
    [ms,Ps,ms_,Ps_,Ds,lh] = SigmaFilter(p0,ys,[],[],[],[],m0,P0);
    [JM,JS] = SigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);
    lhs(k) = lh;
    
    p{5} = as(k);
    [lb,glb] = EM_LB_Ballistic(p,5,ys,JM,JS);

    lbs(k) = lb;
    glbs(k) = glb;
end

figure(1); clf;
subplot(3,1,1); 
plot(as,lhs); grid ON; title('Likelihood');
subplot(3,1,2); 
plot(as,lbs); grid ON; title('Lower bound');
subplot(3,1,3);
plot(as,glbs); grid ON; title('d Lower bound');
figure(2);clf;
plot(as,lhs-lbs);grid on;

% figure(3); clf;
% subplot(3,1,1); 
% plot(as,lhsSR); grid ON; title('LikelihoodSR');
% subplot(3,1,2); 
% plot(as,lbsSR); grid ON; title('Lower bound SR');
% subplot(3,1,3);
% plot(as,glbsSR); grid ON; title('d Lower bound SR');
% figure(4);clf;
% plot(as,lhsSR-lbsSR);grid on;




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





