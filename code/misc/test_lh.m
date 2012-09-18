% test both ways of computing 
% the likelihood and its gradient

%% AR(1) model

N = 400;
a = 0.5;
q = 0.5;
r = 0.01;

x = zeros(1,N);
y = x;

x_ = 5;
for k = 1:N
    x_ = a*x_+sqrt(q)*randn;
    x(k) = x_;
    y(k) = x_ + sqrt(r)*randn;
end

m0 = 0;
P0 = 1;
p0 = {a,q,1,r,0};
[ms,Ps,ms_,Ps_,Ds] = SigmaFilter(p0,y,[],[],[],[],m0,P0);
[JM,JP] = OldSigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);

figure(1); clf;
subplot(2,1,1);
plot(1:N,x,1:N,ms(1,:),1:N,JM(1,:)); grid ON;
subplot(2,1,2);
plot(1:N,squeeze(Ps(1,1,:)),1:N,squeeze(JP(1,1,:))); grid ON;


%% compute likelihood and gradient for A on grid

[eij,wij] = CKFPoints(2);
model = {};
model.eij = eij; model.wij = wij; model.m0 = m0; model.P0 = P0;
model.h = @(x,k,p) x(1,:);
model.f = @(x,k,p) p{1}*x;
model.gradDim = 1;
model.grad = @AR1Grad;

NN = 15;
as = linspace(0.2,0.8,NN);
lhs = zeros(1,NN);
glhs = lhs;
lbs = lhs;
glbs = lhs;
for k=1:NN
    p0 = {as(k),q,1,r,0};
    [ms,Ps,ms_,Ps_,Ds,lh,glh] = SigmaFilter(p0,y,[],[],[],[],m0,P0);
    lhs(k) = lh;
    glhs(k) = glh;
    [JM,JP] = OldSigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);
    [lb,glb] = EM_LB_Sigma(model,p0,y,JM,JP);
    lbs(k) = lb;
    glbs(k) = glb;
end

figure(1); clf;
subplot(2,1,1); 
plot(as,lhs); grid ON;
subplot(2,1,2);
plot(as,glhs,as,glbs); grid ON;

%% compute likelihood and gradient for Q on grid

[eij,wij] = CKFPoints(2);
model = {};
model.eij = eij; model.wij = wij; model.m0 = m0; model.P0 = P0;
model.h = @(x,k,p) x(1,:);
model.f = @(x,k,p) p{1}*x;
model.gradDim = 2;
model.grad = @AR1Grad;

NN = 20;
qs = linspace(0.4,0.6,NN);
lhs = zeros(1,NN);
lbs = lhs;
glbs = lhs;
for k=1:NN
    p0 = {a,qs(k),1,r,0};
    [ms,Ps,ms_,Ps_,Ds,lh,glh] = SigmaFilter(p0,y,[],[],[],[],m0,P0);
    lhs(k) = lh;
    glhs(k) = glh;
    [JM,JP] = OldSigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);
    [lb,glb] = EM_LB_Sigma(model,p0,y,JM,JP);
    lbs(k) = lb;
    glbs(k) = glb(2);
end

figure(1); clf;
subplot(2,1,1); 
plot(qs,lhs); grid ON;
subplot(2,1,2);
plot(qs,glhs,qs,glbs); grid ON;


%% try maximization with BFGS


opt = optimset(@fminunc);
opt.GradObj = 'on';
%opt.Display = 'off';
%opt.OutputFcn = @store_iterations;

lhf = @(yessss) SigmaLH({yessss,q,1,r,0},y,[],[],[],[],m0,P0);
[aa,VAL,EF,OP,GRAD] = fminunc(lhf,0.1,opt)




