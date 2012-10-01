
%% Setup
global dt c
S = load('../data/dataa_villelle.mat','card_data','data_t');
y = S.card_data;
K = S.data_t;
clear S;
offset = 20;
secs = 57;
%secs = 1;
starti = find(K<offset,1,'last'); 
endi = find(K<offset+secs,1,'last');
y = y(starti:endi);
K = K(starti:endi);
K = K-K(1);


ds = round(40*secs/60);
%ds = 15;
% downsample by ds
y = y(1:ds:end)';
K = K(1:ds:end);

dt = K(2)-K(1);
N = length(K);
T = K(end);


% the parameters of this model
r = 0.0001;          % log measurement noise
Nqx = 3;
Nqw = 2;
qx_range = linspace(0.04,0.08,Nqx);
qw_range = linspace(0.2,0.35,Nqw);



[QX,QW] = meshgrid(qx_range,qw_range);


c = 3; % number of harmonics (including the fundamental frequency)
xDim = 2*c+1;
H = [0 repmat([1 0],1,c)];
h = @(x) H*x;
Jh = @(x) H;
f = @(x) sinusoid_f(x);
Jf = @(x) sinusoid_Jf(x);
m0 = [0.3 zeros(1,xDim-1)]';
P0 = eye(xDim);


% adjust K and y to include the zeroth measurement
K = [0 K+dt];
Y = [H*m0 y];


fn = '../data/HarmonicLH_%.0f_%.0f';

SR = sqrt(r);
[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

Z = zeros(Nqw,Nqx);
%tic
for i=1:Nqx
  for j=1:Nqw   
    SQ = sqrt(sinusoid_Q(QW(j,i),QX(j,i),1));
    lh = 0; m = m0; S = P0;
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);

        if k==N+1; break; end; 

      
        [m,S,K,my,CSy] = SigmaKF_Update(m_,S_,Y(:,k+1),h,SR,usig,w);
        %%% CSy and CC are Cholesky decompositions _HERE_ %%%
        Sy = CSy*CSy';
        lh = lh + likelihood(Y(:,k+1)-my,Sy);


    end
    Z(j,i) = lh;
    fprintf('%.2f\n',100*((i-1)*Nqw+j)/(Nqx*Nqw));

  
  end
end
%tm = toc
save(sprintf(fn,Nqx,Nqw));



