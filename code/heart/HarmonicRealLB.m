%% Setup
load('../data/Harmonic_lqw_lqx_10_3751.mat');
Nqx = 100;
Nqw = 100;

% est_em_n(:,10,6)
% ans =
% 
%   -1.168451139502009
%   -2.288206680651371

qx_c = -3.2554;
qw_c = est_em_n(1,10,6);


qx_range = linspace(-3,qx_c,Nqx/2);
dx = qx_range(2)-qx_range(1);
qx_range = [qx_range qx_c+dx:dx:-3.6];
Nqx = numel(qx_range);

qw_range = linspace(-2,qw_c,Nqw/2);
dw = qw_range(2)-qw_range(1);
qw_range = [qw_range qw_c+dw:dw:-0.4];
Nqw = numel(qw_range);



[QX,QW] = meshgrid(qx_range,qw_range);



fn = '../data/HarmonicLB_%.0f_%.0f';


% [usig,w] = CKFPoints(2*c+1);
% w(:,2) = sqrt(w(:,1));
% % compute E-step
% p = p0;
% p(gi) = [qw_c qx_c];
% [lh,~,MM,SS,SQ] = Harmonic_LH(p,Y);
% [MS,SM,DD] = SigmaSmoothSR(MM,SS,f,SQ,usig,w); % D = Smoother Gain
% [I1,I2,I3] = EM_I123_Sigma(f,h,m0,Y,MS,SM,DD);


Z = zeros(Nqw,Nqx);
%tic
jobs = cell(Nqx,1);
for i=1:Nqx
  %jobs{i} = batch(@HarmonicRealInnerLH,1,{QW(:,i),QX(:,i),Nqw,H,Y,SR,usig,w,m0,P0,dt},'PathDependencies','./heart');
  Z(:,i) = HarmonicRealInnerLB(QW(:,i),QX(:,i),Nqw,N,p0,P0,I1,I2,I3);
  i
  %wait(j);
  %diary(j)
  %tmp = getAllOutputArguments(j)
  %Z(:,i) = tmp{1};
  %Z(:,i) = HarmonicRealInnerLH(QW(:,i),QX(:,i),Nqw,H,Y,SR,usig,w,m0,P0,dt);
end
save(sprintf(fn,Nqx,Nqw),'QX','QW','Z');

%%
for i=1:Nqx
  tmp = getAllOutputArguments(jobs{i});
  Z(:,i) = tmp{1};  
end

%%
figure(1); clf;
surf(log(QX),log(QW),Z); hold on;
%plot3([exp(lqx) exp(lqx)],[exp(lqw) exp(lqw)],[min(min(Z)) max(max(Z))]);
xlabel('qx');
ylabel('qw');

%% tm = toc



