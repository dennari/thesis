
%% Setup
global dt H c m0 P0

S = load('../data/dataa_villelle.mat','card_data','data_t');
y = S.card_data;
K = S.data_t;
clear S;
nm = 1;
% take first nm minutes
endi = find(K<nm*60,1,'last');
y = y(1:endi);
K = K(1:endi);


ds = 10;
% downsample by ds
y = y(1:ds:end)';
K = K(1:ds:end);

dt = K(2)-K(1);
N = length(K);
T = K(end);


% plot first 60 seconds
%endi = find(K<60,1,'last');
%plot(K(1:endi),y(1:endi));


% the parameters of this model
lqx = log(0.7);           % log Dynamic model noise spectral density
lqw = log(0.1);           % log angular velocity variance
lr =  log(0.02);          % log measurement noise


c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
f = @(x) sinusoid_f(x);
Q = sinusoid_Q(lqw,repmat(lqx,1,c));              
SQ = chol(Q,'lower');
R = sinusoid_R(lr);
SR = sqrt(R);
m0 = [0.5*2*pi zeros(1,xDim-1)]';
P0 = eye(xDim);


% adjust K and y to include the zeroth measurement
K = [0 K+dt];
Y = [H*m0 y];

%% Plot DRIFTER frequency estimate

S = load('../data/dataa_villelle.mat','card_freq','freq_t');

endi = find(S.freq_t<nm*60,1,'last');
DRIFTER_f = S.card_freq(1:endi);
DRIFTER_ft = S.freq_t(1:endi);
plot(DRIFTER_ft,DRIFTER_f);



%% Test

[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

% add the zeroth measurement


MM = zeros(xDim,N+1); MM(:,1) = m0;
SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
m = m0; lh = 0;  S = P0;

for k=1:(N+1)
    [m_,S_] = SigmaKF_Predict(m,S,f,SQ,[],usig,w);
    
    if k==N+1; break; end; 
    
    y = Y(:,k+1); % notice indexing!!!
    [m,S,~,IM,IS] = SigmaKF_Update(m_,S_,y,h,SR,usig,w);
    
    
    MM(:,k+1) = m;
    SS(:,:,k+1) = S;
    %likelihood(y(:,k+1)-IM,IS)
    lh = lh + likelihood(y-IM,IS);

end
[MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,[],usig,w); 



figure(1); clf;
plot(K,Y,K,H*MM,K,H*MS); grid on;
figure(2); clf;
m = 3;
subplot(m,1,1);
plot(K,sqrt(sum((H*MM-Y).^2,1)),K,sqrt(sum((H*MS-Y).^2,1))); grid on; title('Err');
subplot(m,1,2);
plot(K,MS(1,:)/2/pi,DRIFTER_ft,DRIFTER_f); grid on; title('Base Freq');
subplot(m,1,3);
plot(K,squeeze(abs(SM(1,1,:)))/2/pi); grid on; title('Freq Std');



%% Compute LH and gradient on grid

% parameters are 
% p(1)=lqw,    log angular velocity variance
% p(2)=lr,     log measurement variance
% p(3:3+c-1)   log signal component variances


p0 = [lqw lr repmat(lqx,1,c)];
gi = 1; % which one we're estimating
true = p0(gi);

NN = 25;
lhs = zeros(1,NN); glhs = lhs; glbs = lhs;




rnge = true+log(5);
as = linspace(true-rnge,true+rnge,NN);

p = p0;
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

%save('../../data/simulateHeartR.mat','lhs','glhs','glbs');
%% Try optimization



