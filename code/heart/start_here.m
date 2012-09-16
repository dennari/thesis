
%% Simulate data

N = 500;
T = linspace(0,25,N);            % Time points
dt = T(2)-T(1);

% the parameters of this model
qx = 0.7;              % Dynamic model noise spectral density
qw = 0.4;               % angular velocity variance
r = 0.15;                % measurement noise
x0 = [1;0];             % Initial state

gamma = 0;              % D amping parameter





% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.5*ones(1,cp(1));
 L3 = 2.5*ones(1,N-cp(2));
 x = T((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(T(cp(2))-T(cp(1))))*(x-T(cp(1)))+L1(1);
 fr = [L1 L2 L3];

A = 0.5; 
fr = A*sin(2*pi*0.1*T)+1.5*A;

c = 3; % number of harmonics (including the fundamental frequency)
XD = 2*c+1;
H = [0 repmat([1 0],1,c)];
sim = zeros(2,N,c);
for k = 1:c
    % simulate data
    sim(:,:,k) = simulate_periodic_data(T,[],k*fr,qx,x0,gamma);
end
% measurements with added noise

x1 = sum(sim(1,:,:),3);
y = x1+sqrt(r)*randn(1,numel(T));



h = @(x,k,p) H*x;
f = @(x,k,p) sinusoid_f(x,dt);
Q = sinusoid_Q(qw,[qx qx],dt);              

p0 = {eye(XD),Q,H,r};
m0 = [0.75*2*pi zeros(1,XD-1)]';
P0 = eye(XD);
s = rng;
[ms,Ps,ms_,Ps_,Ds,lh] = SigmaFilter(p0,y,f,h,[],[],m0,P0);
[JM,JP] = SigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);



figure(1); clf;
subplot(3,1,1);
%plot(T,y,'kx',T,sum(sim(1,:,:),3));
filt = sum(ms(2:2:XD,:),1);
smooth = sum(JM(2:2:XD,:),1);
plot(T,x1,T,filt,T,smooth); xlabel('t'); ylabel('x(1)'); grid on;
subplot(3,1,2);
plot(T,fr,T,ms(1,:)/(2*pi),T,JM(1,:)/(2*pi));
title('Signal'); legend('True','Filtered','Smoothed');
subplot(3,1,3);
plot(T,squeeze(Ps(1,1,:)),T,squeeze(JP(6,6,:)));
title('Frequency variance'); legend('Filte','Smoother');


% test square-root filter
p0 = {eye(XD),chol(Q,'lower'),H,sqrt(r)};
rng(s);
[ms,Ss,lh2] = SigmaFilterSR(p0,y,f,h,[],[],m0,P0);
[JM,JS] = SigmaSmootherSR(f,p0{2},ms,Ss,m0,P0);

disp(sqrt((lh-lh2)^2));

figure(2); clf;
subplot(3,1,1);
%plot(T,y,'kx',T,sum(sim(1,:,:),3));
filt = sum(ms(2:2:XD,:),1);
smooth = sum(JM(2:2:XD,:),1);
plot(T,x1,T,filt,T,smooth); xlabel('t'); ylabel('x(1)'); grid on;
subplot(3,1,2);
plot(T,fr,T,ms(1,:)/(2*pi),T,JM(1,:)/(2*pi));
title('Signal'); legend('True','Filtered','Smoothed');
subplot(3,1,3);
plot(T,squeeze(Ss(1,1,:)).^2,T,squeeze(JS(6,6,:)));
title('Frequency variance'); legend('Filter','Smoother');

%plot(T,squeeze(JS(6,6,:))-squeeze(JP(6,6,:)));


%% Compute LH and gradient on grid

% parameters are 
% p{1}=qw, angular velocity variance
% p{2}=qx OR Qx \in 2x2, the signal component variance 
% p{3}=r, measurement variance
% p{4}=m0, prior mean
% p{5}=P0, prior covariance matrix
% p{6}=H, state-to-measurement matrix, fixed
% p{7}=dt, not really a parameter
% set up the starting point
p{1} = qw;
p{2} = 0.1;
p{3} = r;
p{4} = [0.75*2*pi zeros(1,XD-1)]';
p{5} = eye(XD);
p{6} = H;
p{7} = dt;

m0 = p{4}; 
P0 = p{5};
h = @(x,k,p) H*x;
f = @(x,k,p) sinusoid_f(x,dt);

NN = 150;
as = linspace(0.5,0.9,NN);
lhs = zeros(1,NN);
lhsSR = lhs;
glhs = lhs;
lbs = lhs;
lbsSR = lhs;
glbs = lhs;
glbsSR = lhs;
p0 = {eye(XD),[],H,p{3}};
for k=1:NN
    k
    as(k)
    p0{2} = sinusoid_Q(p{1},[as(k) as(k)],dt);
    [ms,Ss,lhSR] = SigmaFilterSR(p0,y,f,h,[],[],m0,P0);
    [JMSR,JSSR] = SigmaSmootherSR(f,p0{2},ms,Ss,m0,P0);
    [ms,Ps,ms_,Ps_,Ds,lh] = SigmaFilter(p0,y,f,h,[],[],m0,P0);
    [JM,JS] = SigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);
    lhs(k) = lh;
    lhsSR(k) = lhSR;
    %glhs(k) = glh;
    
    %
    p{2} = as(k);
    [lb,glb] = em_lb_harmonic(p,2,y,JM,JS);
    [lbSR,glbSR] = em_lb_harmonic(p,2,y,JMSR,JSSR);

    lbs(k) = lb;
    glbs(k) = glb;
    lbsSR(k) = lbSR;
    glbsSR(k) = glbSR;
end

figure(1); clf;
subplot(3,1,1); 
plot(as,lhs); grid ON; title('Likelihood');
subplot(3,1,2); 
plot(as,lbs); grid ON; title('Lower bound');
subplot(3,1,3);
plot(as,glbs); grid ON; title('d Lower bound');
figure(2);clf;
plot(as,lhs-lbs);grid on;

figure(3); clf;
subplot(3,1,1); 
plot(as,lhsSR); grid ON; title('LikelihoodSR');
subplot(3,1,2); 
plot(as,lbsSR); grid ON; title('Lower bound SR');
subplot(3,1,3);
plot(as,glbsSR); grid ON; title('d Lower bound SR');
figure(4);clf;
plot(as,lhsSR-lbsSR);grid on;




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
    