syms t ax ay real;

T = 6.4;
N = 300;
t_ = T/N;
K = (0:N)*t_;
ax_ = (7*t_)^2;
ay_ = 0;
v0 = 40;%150/3.6; % magnitude of the initial velocity
alpha0 = 50/360*2*pi; % initial direction
gy = -9.81;
gx = -3;

L = [t^2/2 t 1 0 0 0;0 0 0 t^2/2 t 1]';
B = diag([ax ay]);
Q = subs(L*B*L',[t ax ay],[t_ ax_ ay_]);
AA = [1 t t^2/2;0 1 t; 0 0 1]; % for one dimension
A = subs(blkdiag(AA,AA),t,t_); % for two dimensions
H = zeros(2,6); H(1,1) = 1; H(2,4) = 1;
vr = 2.7;
R = vr*eye(2);

x0 = [0 v0*cos(alpha0) gx 0 v0*sin(alpha0) gy]';
xs = zeros(6,N+1);
skipY = 12;
ys = zeros(2,N+1);
xs(:,1) = x0;
ys(:,1) = H*x0;


% simulate
x = x0;
for k=1:N
	x = mvnrnd(A*x,Q)';
	xs(:,k+1) = x;
	ys(:,k+1) = NaN;
	if ~mod(k,skipY)
		ys(:,k+1) = mvnrnd(H*x,R)';
	end
end
% the last index where y is still positive
N = find(xs(4,:)>=0,1,'last');
K = (0:(N-1))*t_;
xs = xs(:,1:N);
ys = ys(:,1:N);
yI = ~isnan(ys(1,:));



textwidth = 426.79134/72.27; % latex textwidth in inches
% plot the true locations and the measurements
plt = struct(); ax = subplot(1,1,1);
plt.data = {{xs(1,:) xs(4,:)},{ys(1,yI),ys(2,yI),'x'}};
%plt.xlabel = '$x\,\mathrm{[m]}$';
%plt.ylabel = '$y\,\mathrm{[m]}$';
plt.legend = {'position' 'measurement'};
plt.w = textwidth*0.5;
%plotstruct(ax,plt);
pyplot('../img/ex1_pos_meas.pdf',plt);
plt = struct(); ax = subplot(1,1,1);
plt.data = {{K xs(2,:)},{}};
%plt.xlabel = '$t$';
%plt.ylabel = '$\dot{x}$';
plt.w = textwidth*0.5;
pyplot('../img/ex1_x2.pdf',plt)
plt.data = {{K xs(3,:)},{}};
%plt.ylabel = '$\ddot{x}$';
pyplot('../img/ex1_x3.pdf',plt);
plt.data = {{K xs(5,:)},{}};
%plt.ylabel = '$\dot{y}$';
pyplot('../img/ex1_x5.pdf',plt);
plt.data = {{K xs(6,:)},{}};
%plt.ylabel = '$\ddot{y}$';
pyplot('../img/ex1_x6.pdf',plt);


break

[ms,Ps,ms_,Ps_,Ds] = SigmaFilter({A,Q,H,R},ys,@(x,k,p)p{1}*x,@(x,k,p)p{3}*x,[],[],x0,eye(6));
[mF,PF] = SigmaSmoother(ms,Ps,ms_,Ps_,Ds,x0,eye(6));


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





