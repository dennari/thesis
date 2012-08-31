
%% Simulate resonator

  dt = 0.01;              % Time discrteization
  T = 0:dt:10;            % Time points
  f = 2*ones(1,numel(T)); % Frequency trajectory
  Qc = 0.01;              % Dynamic model noise spectral density
  x0 = [1;0];             % Initial state
  gamma = .9;             % Damping parameter

  % Simulate
  [x] = simulate_periodic_data(T,[],f,Qc,x0,gamma);
  
  % Show
  figure(1); clf
    plot(T,x(1,:))
  
    
%% Load physiological reference data

  % Data file
  refpath = 'drifter_1206200003.txt';
  
  % show the data
  loadReference(refpath);
  
  % or get the data
  [refdata,T] = loadReference(refpath);
    