
%% Setup
global dt H c m0 P0
N = 500;
T = 25;
dt = T/N;

K = (0:N)*dt;

% the parameters of this model
qx = 0.7;              % Dynamic model noise spectral density
qw = 0.1;               % angular velocity variance
r = 0.02;                % measurement noise
R = r;
SR = sqrt(r);

c = 2; % number of harmonics (including the fundamental frequency)
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
MM = zeros(xDim,N+1); MM(:,1) = m0;
SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
dm0 = zeros(xDim,1); dP0 = zeros(xDim);
m = m0; lh = 0;  S = P0;
dm = dm0; dP = dP0;

for k=1:(N+1)
    [m_,S_] = SigmaKF_Predict(m,S,f,SQ,[],usig,w);
    if k==N+1; break; end; 

    [m,S,~,IM,IS] = SigmaKF_Update(m_,S_,y(:,k+1),h,SR,usig,w);
    
    
    MM(:,k+1) = m;
    SS(:,:,k+1) = S;
    %likelihood(y(:,k+1)-IM,IS)
    lh = lh + likelihood(y(:,k+1)-IM,IS);

end
[MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); 



figure(1); clf;
plot(K,H*xs,K,H*MM,K,ys,'kx'); grid on;
figure(2); clf;
m = 3;
subplot(m,1,1);
plot(K,sqrt(sum((MM-xs).^2)),K,sqrt(sum((MS-xs).^2))); grid on; title('Err');
subplot(m,1,2);
plot(K,xs(1,:),K,MM(1,:),K,MS(1,:)); grid on; title('Freq');
subplot(m,1,3);
plot(K,squeeze(abs(SS(1,1,:))),K,squeeze(abs(SM(1,1,:)))); grid on; title('Freq Std');
% figure(2); clf;
% subplot(2,2,1);
% plot(K,xs(1,:)); 
% subplot(2,2,2);
% plot(K,xs(3,:)); 
% subplot(2,2,3);
% plot(K,xs(2,:));
% subplot(2,2,4);
% plot(K,xs(4,:));
% 
% figure(3); clf; n=3;
% subplot(n,1,1);
% plot(K,squeeze(PP(1,1,:)));
% subplot(n,1,2);
% plot(K,squeeze(PP(2,2,:)));
% subplot(n,1,3);
% plot(K,squeeze(PP(3,3,:)));

% figure(1); clf;
% subplot(3,1,1);
% %plot(T,y,'kx',T,sum(sim(1,:,:),3));
% filt = sum(ms(2:2:xDim,:),1);
% smooth = sum(JM(2:2:xDim,:),1);
% plot(T,x1,T,filt,T,smooth); xlabel('t'); ylabel('x(1)'); grid on;
% subplot(3,1,2);
% plot(T,fr,T,ms(1,:)/(2*pi),T,JM(1,:)/(2*pi));
% title('Signal'); legend('True','Filtered','Smoothed');
% subplot(3,1,3);
% plot(T,squeeze(Ps(1,1,:)),T,squeeze(JP(6,6,:)));
% title('Frequency variance'); legend('Filte','Smoother');



%plot(T,squeeze(JS(6,6,:))-squeeze(JP(6,6,:)));


%% Compute LH and gradient on grid

% parameters are 
% p(1)=qw,    angular velocity variance
% p(2)=r,     measurement variance
% p(3:3+c-1)  the signal component variances


p0 = [qw r repmat(qx,1,c)];
p = p0;

NN = 25;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;


gi = 2;

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




%% try harmonic stuff
dt = 1/(9e3); % sample at 20kHz
T = 0:dt:10;            % Time points
N = numel(T);
gamma = 0;
Qc = 0.07;              % Dynamic model noise spectral density
qw = 0.3;               % angular velocity variance
x0 = [1;0];

base = 200/60*200;
fr = base+0.3*base*sin(2*pi*0.2*T);
%figure(1);clf;plot(T,fr); grid ON;
%break
% simulate data
s = rng;
[x1] = simulate_periodic_data(T,[],fr,Qc,x0,gamma);
disp(1)
rng(s);
[x2] = simulate_periodic_data(T,[],2*fr,Qc,x0,gamma);
rng(s);
disp(2)
[x3] = simulate_periodic_data(T,[],3*fr,Qc,x0,gamma);
disp(3)

figure(1);clf;
subplot(4,1,1);
plot(T,x1(1,:)); grid ON;
subplot(4,1,2);
plot(T,x2(1,:)); grid ON;
subplot(4,1,3);
plot(T,x3(1,:)); grid ON;
subplot(4,1,4);
y = 1/3*(x1(1,:)+x2(1,:)+x3(1,:));
plot(T,y); grid ON;
%%
x = [x1(1,:);x2(1,:);x3(1,:)];
coeff = [2 1 1]; coeff = coeff/sum(coeff);
soundsc(x'*coeff',1/dt,16);

%% Load physiological reference data

  % Data file
  %refpath = 'drifter_1206200003.txt';
  
  % show the data
  %loadReference(refpath);
  
  % or get the data
  %[refdata,T] = loadReference(refpath);
    