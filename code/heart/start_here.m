
%% Original model

N = 250;
T = linspace(0,10,N);            % Time points
dt = T(2)-T(1);
A = 0.5;
Qx = 0.1;              % Dynamic model noise spectral density
Qw = 0.1;               % angular velocity variance
x0 = [1;0];             % Initial state
gamma = 0;              % D amping parameter
r = 0.1*2;                % measurement noise


% artificial frequency trajectory
 cp = floor([1 2]*N/3);
 L1 = 0.5*ones(1,cp(1));
 L3 = 2.5*ones(1,N-cp(2));
 x = T((cp(1)+1):cp(2));
 L2 = ((L3(1)-L1(1))/(T(cp(2))-T(cp(1))))*(x-T(cp(1)))+L1(1);
 fr = [L1 L2 L3];

 
 %fr = A*sin(2*pi*0.1*T)+1.5*A;
 
% simulate data
[x] = simulate_periodic_data(T,[],fr,Qx,x0,gamma);
% measurements with added noise
y = x(1,:)+sqrt(r)*randn(1,numel(T));

XD = 3;
h = @(x,k,p) x(1,:);
f = @(x,k,p) sinusoid_A(x(3),dt)*x;
Q = @(m,k,p) sinusoid_Q(m(3),Qx ,Qw,dt);              

m0 = [0.5 0.5 2*pi*0.4]';
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
    