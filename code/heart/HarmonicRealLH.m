%% Setup
load('../data/Harmonic_lqw_lqx_10_3751.mat');
Nqx = 50;
Nqw = 50;
qx_range = linspace(exp(-5.5),exp(-1.5),Nqx);
qw_range = linspace(exp(-2.2),exp(-0.05),Nqw);

[QX,QW] = meshgrid(qx_range,qw_range);



fn = '../data/HarmonicLH_%.0f_%.0f';

SR = sqrt(sinusoid_R(lr));
[usig,w] = CKFPoints(xDim);
w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth

Z = zeros(Nqw,Nqx);
%tic
jobs = cell(Nqx,1);
for i=1:Nqx
  jobs{i} = batch(@HarmonicRealInnerLH,1,{QW(:,i),QX(:,i),Nqw,H,Y,SR,usig,w,m0,P0,dt},'PathDependencies','./heart');
  %wait(j);
  %diary(j)
  %tmp = getAllOutputArguments(j)
  %Z(:,i) = tmp{1};
  %Z(:,i) = HarmonicRealInnerLH(QW(:,i),QX(:,i),Nqw,H,Y,SR,usig,w,m0,P0,dt);
end

%%
for i=1:Nqx
  tmp = getAllOutputArguments(jobs{i});
  Z(:,i) = tmp{1};  
end

%%

surf(QX,QW,Z); hold on;
%plot3([exp(lqx) exp(lqx)],[exp(lqw) exp(lqw)],[min(min(Z)) max(max(Z))]);
xlabel('qx');
ylabel('qw');

%% tm = toc
save(sprintf(fn,Nqx,Nqw),'QX','QW','Z');


