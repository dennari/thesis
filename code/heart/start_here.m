
%% Setup
global dt c m0 P0 h f Jh Jf
N = 500;
T = 15;
dt = T/N;
K = (0:N)*dt;

% the parameters of this model
lqx = log(0.6);    % log(sqrt) Dynamic model noise spectral density
lqw = log(0.5);   % log(sqrt) angular velocity noise variance
lr =  log(0.05);   % log(sqrt) measurement noise



c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x);
Jf = @(x) sinusoid_Jf(x);
Q = sinusoid_Q(lqw,lqx);              

if Q == diag(diag(Q))
  SQ = sqrt(Q);
else
  SQ = chol(Q,'lower');
end
R = sinusoid_R(lr);
SR = chol(R,'lower');
m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


%% Simulate


% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.8*ones(1,cp(1));
 L3 = 1.0*ones(1,N-cp(2)+1);
 x = K((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(K(cp(2))-K(cp(1))))*(x-K(cp(1)))+L1(1);
 fr = 2*pi*[L1 L2 L3];
%fr = 2*pi*0.1*ones(1,N+1);

%A = 0.5; 
%fr = A*sin(2*pi*0.1*K)+1.5*A;
%fr = ones(size(K))*2*pi*0.8;

[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
xs = zeros(xDim,N+1);
ys = zeros(1,N+1);

x = m0;
x(1) = 0.5*2*pi;
xs(:,1) = x;
ys(:,1) = h(x);

for k=2:N+1
	x = mvnrnd(f(x),Q)';
  %x(1) = fr(k);
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
    [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);
    if k==N+1; break; end; 

    [m,S,~,my,Sy] = SigmaKF_Update(m_,S_,y(:,k+1),h,SR,usig,w);
    
    
    MM(:,k+1) = m;
    SS(:,:,k+1) = S;
    %likelihood(y(:,k+1)-IM,IS)
    lh = lh + likelihood(y(:,k+1)-my,Sy*Sy');

end
[MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); 



figure(1); clf;
plot(K,H*xs,K,H*MM,K,H*MS,K,ys,'kx'); grid on;
figure(2); clf;
m = 3;
subplot(m,1,1);
plot(K,sqrt(sum((MM-xs).^2)),K,sqrt(sum((MS-xs).^2))); grid on; title('Err');
subplot(m,1,2);
plot(K,xs(1,:)/(2*pi),K,MM(1,:)/(2*pi),K,MS(1,:)/(2*pi)); grid on; title('Freq');
subplot(m,1,3);
plot(K,squeeze(abs(SS(1,1,:))),K,squeeze(abs(SM(1,1,:)))); grid on; title('Freq Std');



%% Compute LH and gradient on grid

% parameters are 
% p(1)=lqw,    log sqrt angular velocity variance
% p(2)=lr,     log sqrt measurement variance
% p(3:3+c-1)   log sqrt signal component variances


p0 = [lqw lr lqx];
gi = 3; % which one we're estimating
true = p0(gi);

NN = 25;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;




start = true + log(0.8);
endd = true - log(0.8);
as = linspace(start,endd,NN);
%as = log(linspace(0.06,0.09,NN));

p = p0;
for j=1:NN
    j
    p(gi) = as(j);
    
    [lh,glh,MM,SS,SQ] = Harmonic_LH(p,ys,gi);
    glh
    lhs(j) = lh;
    glhs(j) = glh;
    
    [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); % D = Smoother Gain
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,ys,MS,SM,DD);
    glbs(j) = EM_LB_Harmonic(p,MS(:,1),gi,N,I1,I2,I3);
end




%load('../data/HarmonicTesting.mat');


n = 4; m= 1; %true = exp(true);
figure(1); clf; eas = exp(as); etr = exp(true);
lhax = subplot(n,m,1);
plot(eas,lhs'); grid on; title('likelihood'); xl = xlim; 
%hold on;plot([etr etr],ylim,'-r');  xl = xlim;
dnumax = subplot(n,m,2);
plot(eas(2:end),diff(lhs)./diff(as)); grid on; title('dNum'); xlim(xl);
%hold on;plot([etr etr],ylim,'-r'); 
dsensax = subplot(n,m,3);
plot(eas,glhs); grid on; title('dSens'); xlim(xl);
%hold on;plot([etr etr],ylim,'-r'); 
demax = subplot(n,m,4);
plot(eas,glbs);grid on; title('dEM'); xlim(xl); 


figure(2); clf; 
plot(eas(2:end),diff(lhs)./diff(as),eas,glhs,eas,glbs); grid on;


% figure(3); clf; plot(eas,glbs);grid on; title('dEM')
% %hold on;plot([etr etr],ylim,'-r'); 
% drawnow;
% ylim(dsensax,ylim(dnumax)); ylim(demax,ylim(dnumax));


%% Joint LH






    
