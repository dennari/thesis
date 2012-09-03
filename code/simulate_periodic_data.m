function [x] = simulate_periodic_data(N,dt,f,Qc,x0,gamma,Ff)
%% simulate_and_estimate - Estimate results for one draw 
% 
% Syntax:
%   simulate_periodic_data(N,dt,f,Qc,x0,g)
%
% In:
%   N  - Number of steps to simulate
%   dt - Time discrteization step length (in seconds)
%   f  - A vector of frequencies in Hz
%   Qc - Spectral density of the dynamic noise term
%   x0 - Initial state (default: random)
%   g  - Damping constant (default: 0)
%
% Out:
%   x  - 2xN-vector of simulated states
%
% Description:
%   Simulate a stochastic oscillator with a given frequency trajectory. 
%   The model is set up as a continuous-time state space model, or a 
%   stochastic differential equation. Runnig this code requires the 
%   'lti_disc' function that is available i.e. in the DRIFTER toolbox.
%
% See also:
%   lti_disc 
%   
%   Copyright 2012 Arno Solin, last updated June 18,2012

%% Simulate a stochastic oscillator
  
  % Alternate inputs; instead of N a vector of time steps is given
  if numel(N)>1 && isempty(dt)
    dt = N(2)-N(1);
    N  = numel(N);
  end

  % Allocate space for results
  x = zeros(2,N);
  
  % Initial state
  if nargin < 5 || isempty(x0)
    x(:,1) = randn(2,1);
  else
    x(:,1) = x0;  
  end

  % Damping factor
  if nargin < 6 || isempty(gamma)
    gamma = 0;
  end
  if nargin < 7
    Ff =@(f,g) [ 0         2*pi*f; 
                -2*pi*f   -g];
  end
  
  % The rest of the states
  L = [0;1];
  
  for k=2:N
      
    % Dynamic model

    F = Ff(f(k),gamma);

    
    % Discretize
    [A,Q] = splti_disc(F,L,Qc,dt);
    
    % Determine if stochastic oscillator or deterministic
    if Qc<eps
      x(:,k) = A*x(:,k-1);  
    else
      x(:,k) = A*x(:,k-1) + chol(Q)'*randn(2,1);
    end
      
  end
  