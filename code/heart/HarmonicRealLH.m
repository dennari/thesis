
%% Setup
load('../data/Harmonic_lqw_lqx_10_3751.mat');
Nqx = 50;
Nqw = 50;
qx_range = linspace(exp(-5.5),exp(-1.5),Nqx);
qw_range = linspace(exp(-2.2),exp(-0.05),Nqw);



[QX,QW] = meshgrid(qx_range,qw_range);
%break

% c = 3; % number of harmonics (including the fundamental frequency)
% xDim = 2*c+1;
% H = [0 repmat([1 0],1,c)];
% h = @(x) H*x;
% f = @(x) sinusoid_f(x);
% m0 = [0.3 zeros(1,xDim-1)]';
% P0 = eye(xDim);
% 
% 
% % adjust K and y to include the zeroth measurement
% K = [0 K+dt];
% Y = [H*m0 y];
if matlabpool('size') == 0 % checking to see if my pool is already open
    matlabpool open 20
end

fn = '../data/HarmonicLH_%.0f_%.0f';

SR = sqrt(sinusoid_R(lr));
[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
Qf = @(qw,qx) diag([qw^2, 0, qx^2, 0, qx^2, 0, qx^2]);
f = @(x) sinusoid_f(x,dt);
Z = zeros(Nqw,Nqx);
%tic
parfor i=1:Nqx
  tmp = zeros(Nqw,1);
  for j=1:Nqw   
    SQ = sqrt(Qf(QW(j,i),QX(j,i)));
    lh = 0; m = m0; S = P0;
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);

        if k==N+1; break; end; 

      
        [m,S,K,my,CSy] = SigmaKF_Update(m_,S_,Y(:,k+1),h,SR,usig,w);
        %%% CSy and CC are Cholesky decompositions _HERE_ %%%
        Sy = CSy*CSy';
        lh = lh + likelihood(Y(:,k+1)-my,Sy);


    end
    tmp(j) = lh;
    %fprintf('%.2f\n',100*((i-1)*Nqw+j)/(Nqx*Nqw));

  
  end
  i
  Z(:,i) = tmp;
end
%tm = toc
save(sprintf(fn,Nqx,Nqw));




