% test both ways of computing 
% the likelihood and its gradient

%% AR(2) model

N = 400;
a1 = [0.5 0.5]; a1 = a1/(sum(a1)*(1.02));
A = [a1; 
     1   0];
q1 = 0.3;
q2 = q1/100;
Q = [q1 0;0 q2];
%R = diag([q1*0.7 q1*5 q1*5]);
R = q1*0.7;
%H = [1 0;0 1;1 1];
H = [1 0];
x = zeros(2,N);
y = zeros(1,N);

[usig,w] = CKFPoints(2);

xDim = size(Q,1);
yDim = size(R,1);
m0 = zeros(2,1);
P0 = eye(2);
x_ = mvnrnd(m0,P0)';
x(:,1) = x_;
for k = 1:N
    x_ = A*x_+mvnrnd(zeros(xDim,1),Q)';
    x(:,k+1) = x_;
    y(:,k) = H*x_ + mvnrnd(zeros(yDim,1),R)';
end

f = @(x) A*x;
h = @(x) H*x;

MM1 = zeros(2,N+1);
PP1 = zeros(2,2,N+1);
MM1_ = zeros(2,N+1);
PP1_ = zeros(2,2,N+1);
MM2 = zeros(2,N+1);
PP2 = zeros(2,2,N+1);
MM2_ = zeros(2,N+1);
PP2_ = zeros(2,2,N+1);
CC = zeros(2,2,N+1);

m1 = m0; MM1(:,1) = m0; 
P1 = P0; PP1(:,:,1) = P0; 
m2 = m0; MM2(:,1) = m0;
P2 = P0; PP2(:,:,1) = P0;
Qsr = chol(Q,'lower');
Rsr = chol(R,'lower');
for k=1:(N+1) % playing with the indices to get them to match
    [m1,P1] = kf_predict(m1,P1,A,Q);
    %[m1,P1] = kfsr_predict(m1,P1,A,Qsr);
    MM1_(:,k) = m1;
    PP1_(:,:,k) = P1;
    
    [m2,P2,C] = SigmaKF_Predict(m2,P2,f,Q,usig,w);
    MM2_(:,k) = m2;
    PP2_(:,:,k) = P2;
    CC(:,:,k) = C;
    %rmse(m2,m1)
    if k==N+1; break; end; % ends to prediction after last measurement
    
    [m1,P1] = kf_update(m1,P1,y(:,k),H,R);
    %[m1,P1] = kfsr_update(m1,P1,y(:,k),H,Rsr);
    MM1(:,k+1) = m1;
    PP1(:,:,k+1) = P1;

    [m2,P2] = SigmaKF_Update(m2,P2,y(:,k),h,R,usig,w);
    MM2(:,k+1) = m2;
    PP2(:,:,k+1) = P2;
    %rmse(m2,m1)

end
[MS1,PS1,D,DD] = rts_smooth(MM1,PP1,A,Q);
%[MS1,PS1,DD] = rtssr_smooth(MM1,PP1,A,Q);
[MS2,PS2,DD2] = SigmaSmoother(MM2,PP2,MM2_,PP2_,CC);

DDv = reshape(DD,numel(DD),1);
DD2v = reshape(DD2,numel(DD2),1);
fprintf(1,'RTS:  %.15f\n',rmse(x,MS1));
fprintf(1,'CRTS: %.15f\n',rmse(x,MS2));
fprintf(1,'CROSS:  %.15f\n',rmse(DDv,DD2v));
%fprintf(1,'SR-CRTS: %.5f\n',rmse(x,JMSR(3:4,:)));
%% Plot

% Cubature Kalman filter
figure(1); clf;
subplot(3,1,1);
plot(1:N+1,x(1,:),1:N+1,MM1(1,:),1:N+1,MS1(1,:)); grid ON;
subplot(3,1,2);
plot(1:N+1,x(2,:),1:N+1,MM1(2,:),1:N+1,MS1(2,:)); grid ON;
subplot(3,1,3);
plot(1:N+1,squeeze(PP1(1,1,:)),1:N+1,squeeze(PS1(1,1,:))); grid ON;


%% compute likelihood and gradient for A on grid


% reset
A = [a1; 
     1   0];
Q = [q1 0;0 q2];
%R = diag([q1*0.7 q1*5 q1*5]);
R = q1*0.7;
%H = [1 0;0 1;1 1];
H = [1 0];



NN = 15;
as = linspace(0.1,0.5,NN); % for q1
%as = linspace(0.05,0.3,NN); % for r

lhs = zeros(1,NN);
glhs = lhs;
glbs = lhs;
vi = 3;
MM = zeros(2,N+1); MM_ = MM;
MM(:,1) = m0;
PP = zeros(2,2,N+1); PP_ = PP;
PP(:,:,1) = P0;
CC = PP;
dm0 = m0;
dP0 = zeros(size(P0));

disp('-----------------------------------------------');
for j=1:NN
    j
    Q = [as(j) 0;0 q2];
    %R = as(j);
    m = m0; P = P0; lh = 0;glh = 0;dm = dm0; dP = dP0;
    for k=1:(N+1)
    
        %[m,P,C] = SigmaKF_Predict(m,P,f,Q,usig,w);
        [m,P] = kf_predict(m,P,A,Q);
        MM_(:,k) = m;
        PP_(:,:,k) = P;
        %CC(:,:,k) = C;
    
        if k==N+1; break; end; 

        %[m,P,C,S,d] = SigmaKF_Update(m,P,y(:,k),h,R,usig,w);
        [m,P,IM,IS] = kf_update(m,P,y(:,k),H,R);
        MM(:,k+1) = m;
        PP(:,:,k+1) = P;
        lh = lh + likelihood(y(:,k)-IM,IS);
        [dm,dP,GLH] = AR2_sensitivity(dm,dP,A,H,vi,y(:,k)-IM,IS,PP_(:,:,k)*H',...
                      MM(:,k),PP(:,:,k));
        glh = glh + GLH; 
    end
    glhs(j) = glh;
    lhs(j) = lh;
    
    
    
    %glhs(k) = glh;
    %[MS,PS,DD] = SigmaSmoother(MM,PP,MM_,PP_,CC);
    [MS,PS,D,DD] = rts_smooth(MM,PP,A,Q);
    
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,y,MS,PS,DD);
       
    %[I1,I2,I3] = EM_I123(A,H,m0,y,MS,PS,DD);
    glbs(j) = EM_LB_AR2(vi,Q,R,P0,m0,MS(:,1),N,I1,I2,I3);
    %lbs(2,j) =   EM_LB(P0,Q,R,N,I1,I2,I3);
    %break;
    %glbs(k) = glb;
end
%break;
n = 2;
figure(1); clf;
subplot(n,1,1);
plot(as,lhs); grid on;
subplot(n,1,2); 
plot(as,glhs,as,glbs); grid ON;
%subplot(n,1,3); 
%plot(as,lhs-lbs(2,:)); grid ON;
%subplot(4,1,4); 
%plot(as,lbs(2,:)-lbs(1,:)); grid ON;

%figure(2); clf;
%plot(as,lbs(2,:)); grid ON;
% 
% figure(2); clf;
% plot(as,glbs); grid ON;

%% compute likelihood and gradient for Q on grid

[eij,wij] = CKFPoints(2);
model = {};
model.eij = eij; model.wij = wij; model.m0 = m0; model.P0 = P0;
model.h = @(x,k,p) x(1,:);
model.f = @(x,k,p) p{1}*x;
model.gradDim = 2;
model.grad = @AR1Grad;

NN = 100;
qs = linspace(0.2,0.8,NN);
lbs = lhs;
glbs = lhs;
for k=1:NN
    p0 = {a,qs(k),1,R,0};
    [ms,Ps,ms_,Ps_,Cs,lh,glh] = SigmaFilter(p0,y,[],[],[],[],m0,P0);
    [JM,JP] = SigmaSmoother(ms,Ps,ms_,Ps_,Cs,m0,P0);
    [lb,glb] = EM_LB_Sigma(model,p0,y,JM,JP);
    lbs(k) = lb;
    glbs(k) = glb(2);
end

figure(1); clf;
subplot(2,1,1); 
plot(qs,lbs); grid ON;
subplot(2,1,2);
plot(qs,glbs); grid ON;


%% try maximization with BFGS


opt = optimset(@fminunc);
opt.GradObj = 'on';
%opt.Display = 'off';
%opt.OutputFcn = @store_iterations;

lhf = @(yessss) SigmaLH({yessss,q,1,R,0},y,[],[],[],[],m0,P0);
[aa,VAL,EF,OP,GRAD] = fminunc(lhf,0.1,opt)




