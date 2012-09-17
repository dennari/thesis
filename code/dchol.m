function [D,L] = dchol(P,dP)
%DCHOL  Derivative of Cholesky decomposition
%
% Syntax:
%   [D,L] = dchol(P,dP)
%
% In:
%   P  - Symmetric positive definite matrix
%   dP - Derivative dP/dt of the matrix
%
% Out:
%   D - Derivative dL/dt of Choleky factor of P
%   L - Cholesky factor of P
%
% Description:
%
%   Calculate derivative of Cholesky factor L of P
%   such that P = L L'. That is, given P and dP/dt
%   calculate dL/dt.
%
% Example:
%   dt=1e-6;
%   L1=chol(P+dt*dP)';
%   L0=chol(P)';
%   (L1-L0)/dt   % Numerical derivative
%
%   dchol(P,dP)  % Exact derivative

% Copyright (C) 2004 Simo Särkkä

  L = zeros(size(P));
  D = zeros(size(P));

  for i=1:size(L,1)
    L(i,i) = sqrt(P(i,i) - sum(L(i,1:i-1).^2));
    D(i,i) = (dP(i,i) - 2 * sum(D(i,1:i-1).*L(i,1:i-1)))/L(i,i)/2;
    for j=i+1:size(L,1)
      L(j,i) = (P(j,i) - sum(L(j,1:i-1).*L(i,1:i-1))) / L(i,i);
      D(j,i) = (dP(j,i) - sum(D(j,1:i-1).*L(i,1:i-1) ...
             + L(j,1:i-1).*D(i,1:i-1))) / L(i,i) - (D(i,i)/L(i,i))*L(j,i);
    end
  end

