
%% Setup
global dt H c m0 P0
N = 2500;
T = 25;
dt = T/N;

K = (0:N)*dt;

% the parameters of this model
qx = 0.7;              % Dynamic model noise spectral density
qw = 0.1;               % angular velocity variance
r = 0.02;                % measurement noise
R = r;
SR = sqrt(r);

c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];

h = @(x) H*x;
f = @(x) sinusoid_f(x);
Q = sinusoid_Q(qw,repmat(qx,1,c));              
SQ = chol(Q,'lower');

m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.5*ones(1,cp(1));
 L3 = 2.5*ones(1,N-cp(2)+1);
 x = K((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(K(cp(2))-K(cp(1))))*(x-K(cp(1)))+L1(1);
 fr = 2*pi*[L1 L2 L3];

%A = 0.5; 
%fr = A*sin(2*pi*0.1*K)+1.5*A;
%fr = ones(size(K))*2*pi*0.8;

[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
xs = zeros(xDim,N+1);
ys = zeros(1,N+1);

m0(1) = fr(1); 
x0 = m0;
x = x0;
xs(:,1) = x0;
ys(:,1) = H*x0;

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  x(1) = fr(k);
  xs(:,k) = x;
	ys(:,k) = mvnrnd(h(x),R)';
  
end
%u = zeros(size(x));
y = ys;

%plot(T,squeeze(JS(6,6,:))-squeeze(JP(6,6,:)));


%% Compute LH and gradient on grid

% parameters are 
% p(1)=qw,    angular velocity variance
% p(2)=r,     measurement variance
% p(3:3+c-1)  the signal component variances


p0 = [qw r repmat(qx,1,c)];
p = p0;

NN = 200;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;


gi = 1;

as = linspace(0.5*p0(gi),1.5*p0(gi),NN);


for j=1:NN
    j
    p(gi) = as(j);
    
    [lh,glh,MM,SS] = Harmonic_LH(p,ys,gi);
    lhs(j) = lh;
    glhs(j) = glh;
    [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); % D = Smoother Gain
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,ys,MS,SM,DD);
    glbs(j) = EM_LB_Harmonic(p,MS(:,1),gi,N,I1,I2,I3);
end

n = 3; m= 1;
figure(1); clf;
subplot(n,m,1);
plot(as,lhs'); grid on;
subplot(n,m,2);
plot(as,glhs,as,glbs); grid on;
subplot(n,m,3);
plot(as,sqrt((glbs-glhs).^2)); grid on;

save('../../data/simulateHeartR.mat','lhs','glhs','glbs');

    
