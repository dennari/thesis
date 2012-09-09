% test both ways of computing 
% the likelihood and its gradient

%% AR(2) model

N = 400;
a1 = [0.5 0.5]; a1 = a1/(sum(a1)*(1+1e-5));
A = [a1; 
     1   0];
q1 = 0.3;
Q = [q1 0;0 0];
R = [q1*0.7 0;0 q1*5];
H = [1 0;0 1];
x = zeros(2,N);
y = zeros(2,N);

[usig,w] = CKFPoints(2);


m0 = zeros(2,1);
P0 = eye(2);
x_ = mvnrnd(m0,P0)';
x(:,1) = x_;
for k = 1:N
    x_ = A*x_+mvnrnd(zeros(2,1),Q)';
    x(:,k+1) = x_;
    y(:,k) = H*x_ + mvnrnd(zeros(2,1),R)';
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
for k=1:(N+1) % playing with the indices to get them to match
    [m1,P1] = kf_predict(m1,P1,A,Q);
    MM1_(:,k) = m1;
    PP1_(:,:,k) = P1;
    
    [m2,P2,C] = SigmaKF_Predict(m2,P2,f,Q,usig,w);
    MM2_(:,k) = m2;
    PP2_(:,:,k) = P2;
    CC(:,:,k) = C;
    
    if k==N+1; break; end; % ends to prediction after last measurement
    
    [m1,P1] = kf_update(m1,P1,y(:,k),H,R);
    MM1(:,k+1) = m1;
    PP1(:,:,k+1) = P1;

    [m2,P2] = SigmaKF_Update(m2,P2,y(:,k),h,R,usig,w);
    MM2(:,k+1) = m2;
    PP2(:,:,k+1) = P2;

end
[MS1,PS1] = rts_smooth(MM1,PP1,A,Q);
[MS2,PS2] = SigmaSmoother(MM2,PP2,MM2_,PP2_,CC);


fprintf(1,'RTS:  %.12f\n',rmse(x,MS1));
fprintf(1,'CRTS: %.12f\n',rmse(x,MS2));
%fprintf(1,'SR-CRTS: %.5f\n',rmse(x,JMSR(3:4,:)));
%% Plot

% Cubature Kalman filter
figure(1); clf;
subplot(3,1,1);
plot(1:N,x(1,:),1:N,ms(1,:),1:N,JM(3,:)); grid ON;
subplot(3,1,2);
plot(1:N,x(2,:),1:N,ms(2,:),1:N,JM(4,:)); grid ON;
subplot(3,1,3);
plot(1:N,squeeze(Ps(1,1,:)),1:N,squeeze(JP(3,1,:))); grid ON;

% SR Cubature Kalman filter
figure(2); clf;
subplot(3,1,1);
plot(1:N,x(1,:),1:N,ms(1,:),1:N,JMSR(3,:)); grid ON;
subplot(3,1,2);
plot(1:N,x(2,:),1:N,ms(2,:),1:N,JMSR(4,:)); grid ON;
subplot(3,1,3);
plot(1:N,squeeze(Ps(1,1,:)),1:N,squeeze(JP(3,1,:)),1:N,squeeze(Ss(1,1,:)).^2,1:N,squeeze(JS(3,1,:))); grid ON;


%% compute likelihood and gradient for A on grid

NN = 150;
as = linspace(0.1,0.5,NN);
lhs = zeros(1,NN);
glhs = lhs;
lbs = lhs;
glbs = lhs;

MM = zeros(2,N+1); MM_ = MM;
PP = zeros(2,2,N+1); PP_ = PP;
CC = PP;

for j=1:NN
    Q = [as(j) 0;0 0];
    m = m0; P = P0;
    lh = 0;
    for k=1:(N+1)
    
        [m,P,C] = SigmaKF_Predict(m,P,f,Q,usig,w);
        MM_(:,k) = m;
        PP_(:,:,k) = P;
        CC(:,:,k) = C;
    
        if k==N+1; break; end; 

        [m,P,S,d] = SigmaKF_Update(m,P,y(:,k),h,R,usig,w);
        MM(:,k+1) = m;
        PP(:,:,k+1) = P;
        lh = lh + likelihood(d,S);
    end
    lhs(j) = lh;
    %glhs(k) = glh;
    [MS,PS,DD] = SigmaSmoother(MM,PP,MM_,PP_,CC);
    [I1,I2,I3] = EM_I123_Sigma(f,h,m0,Y,MS,PS,DD);
    lbs(k) =     EM_LB(P0,Q,R,N,I1,I2,I3);
    %glbs(k) = glb;
end

figure(1); clf;
plot(as,lhs); grid on;
% subplot(2,1,1); 
% plot(as,lhs,as,lbs+(lhs(2)-lbs(2))); grid ON;
% subplot(2,1,2);
% plot(as,glhs,as,glbs); grid ON;

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




