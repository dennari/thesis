
%% Simulate resonator

  dt = 0.02;              % Time discrteization
  T = 0:dt:100;            % Time points
  f = 1.2*ones(1,numel(T)); % Frequency trajectory
  Qc = 0.01;              % Dynamic model noise spectral density
  Qw = 0.1;
  x0 = [1;0];             % Initial state
  gamma = 0;             % D amping parameter



% simulate the frequency from AR(1) model
%xf = 5;
%f(1) = xf;
%for k=2:numel(T)
%  xf = 0.999*xf+0.01*2*randn;
%  f(k) = xf;
%end

A = 2;
f = A*sin(2*pi*0.1*T)+1.2*A;
% k = 0.7;
% b = 0.7;
% L1 = k*T(1:500)+b;
% L2 = -k*T(501:end)+k*(T(500)+T(501))+b;
% f = [L1 L2];

%figure(1); clf;
%plot(T,f);
%pause
% continuous time transition matrix
% F = @(f,g) [0 1; -(2*pi*f)^2 -g];

[x] = simulate_periodic_data(T,[],f,Qc,x0,gamma);
r = 0.1;
% measurements with added noise
y = x(1,:)+sqrt(r)*randn(1,numel(T));  
  % Show
  figure(1); clf; 
  subplot(3,1,1);
  plot(T,x(1,:),T,y,'xk')
  %subplot(3,1,2);
  %plot(T,x(2,:))

XD = 3;
m0 = [0 0 3]';
P0 = eye(XD);
% fake linear parameters
p0 = {eye(XD),eye(XD),[1 0 0],r};

F = @(w) [ 0   2*pi*w      0; 
          -2*pi*w  -gamma  0;
           0   0      0];
L = [0 1 0; 0 0 1]';
disc = @(k,m,P) splti_disc(F(m(3)),L,diag([Qc Qw]),dt);

[ms,Ps,ms_,Ps_,Ds] = SigmaFilterPeriodic(p0,y,[],@(x,k,p)p{3}*x,[],[],m0,P0,disc);
subplot(3,1,2);
plot(T,x(1,:),T,ms(1,:))
subplot(3,1,3);
plot(T,f,T,ms(3,:))

%% Try model from Kim et al. 2010

dt = 0.02;              % Time discrteization
T = 0:dt:10;            % Time points
A = 2;
fr = A*sin(2*pi*0.1*T)+1.2*A;
Qc = 0.01;              % Dynamic model noise spectral density
Qw = 0.1;
x0 = [1;0];             % Initial state
gamma = 0;              % D amping parameter
r = 0.1;                % measurement noise
N = numel(T);



% simulate data
[x] = simulate_periodic_data(T,[],fr,Qc,x0,gamma);
% measurements with added noise
y = x(1,:)+sqrt(r)*randn(1,numel(T));

XD = 3;
h = @(x,k,p) x(2,:).*cos(2*pi*dt*k*x(1,:))+x(3,:).*sin(2*pi*dt*k*x(1,:));
f = @(x,k,p) x;
% frequency, cos-component, sin-component
m0 = [1 0 0]';
P0 = eye(XD);
Q = diag([0.1 0.1 0.01]);
p0 = {[],Q,[],r};
[ms,Ps,ms_,Ps_,Ds] = SigmaFilter(p0,y,f,h,[],[],m0,P0);
  
figure(1); clf;
subplot(2,1,1);
plot(T,sqrt(ms(2,:).^2+ms(3,:).^2),T,x(1,:))
subplot(2,1,2);
plot(T,fr,T,ms(1,:))

%% Original model

N = 250;
T = linspace(0,25,N);            % Time points
dt = T(2)-T(1);
A = 0.5;
Qc = 0.07;              % Dynamic model noise spectral density
Qw = 0.3;               % angular velocity variance
x0 = [1;0];             % Initial state
gamma = 0;              % D amping parameter
r = 0.1*2;                % measurement noise


% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.5*ones(1,cp(1));
 L3 = 2.5*ones(1,N-cp(2));
 x = T((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(T(cp(2))-T(cp(1))))*(x-T(cp(1)))+L1(1);
% fr = [L1 L2 L3];
fr = A*sin(2*pi*0.1*T)+1.5*A;
 
% simulate data
[x] = simulate_periodic_data(T,[],fr,Qc,x0,gamma);
% measurements with added noise
y = x(1,:)+sqrt(r)*randn(1,numel(T));

XD = 3;
h = @(x,k,p) x(1,:);
f = @(x,k,p) sinusoid_f(x,k,p,dt);
Q = @(m,k,p) sinusoid_Q(m,k,p,Qc ,Qw,dt);              

m0 = [0.5 0.5 2*pi*0.5]';
P0 = eye(XD)*0.1;

p0 = {eye(XD),Q,[1 0 0],r};
[ms,Ps,ms_,Ps_,Ds] = SigmaFilter(p0,y,f,h,[],[],m0,P0);
[JM,JP] = SigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0);


figure(1); clf;
subplot(2,1,1);
%plot(T,y,'kx',T,x(1,:),T,ms(1,:))
plot(T,x(1,:),T,ms(1,:),T,JM(1,:)); xlabel('t'); ylabel('x(1)');
title('Signal'); legend('True','Filtered','Smoothed');
subplot(2,1,2);
plot(T,fr,T,abs(ms(3,:))/(2*pi),T,abs(JM(3,:))/(2*pi));xlabel('t'); ylabel('x(3)');
title('Frequency'); legend('True','Filtered','Smoothed');

%% try harmonic stuff
dt = 1/(9e3); % sample at 20kHz
T = 0:dt:10;            % Time points
N = numel(T);
gamma = 0;
Qc = 0.07;              % Dynamic model noise spectral density
Qw = 0.3;               % angular velocity variance
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
    