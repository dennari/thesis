%RTS_SMOOTH  Rauch-Tung-Striebel smoother
%
% Syntax:
%   [M,P,S] = RTS_SMOOTH(M,P,M_,P_)
%
% In:
%   M - NxK matrix of K mean estimates from Kalman filter
%   P - NxNxK matrix of K state covariances from Kalman Filter
%   M_- NxK matrix of K prediction mean estimates from Kalman filter
%   P_- NxNxK matrix of K prediction state covariances from Kalman Filter
%   A - NxN state transition matrix or NxNxK matrix of K state
%       transition matrices for each step.
%
% Out:
%   M - Smoothed state mean sequence
%   P - Smoothed state covariance sequence
%   D - Smoother gain sequence
%   
% Description:
%   Rauch-Tung-Striebel smoother algorithm. Calculate "smoothed"
%   sequence from given Kalman filter output sequence
%   by conditioning all steps to all measurements.
%
% Example:
%   m = m0;
%   P = P0;
%   MM = zeros(size(m,1),size(Y,2));
%   PP = zeros(size(m,1),size(m,1),size(Y,2));
%   for k=1:size(Y,2)
%     [m,P] = kf_predict(m,P,A,Q);
%     [m,P] = kf_update(m,P,Y(:,k),H,R);
%     MM(:,k) = m;
%     PP(:,:,k) = P;
%   end
%   [SM,SP] = rts_smooth(MM,PP,A,Q);
%
% See also:
%   KF_PREDICT, KF_UPDATE

% Copyright (C) 2003-2006 Simo S�rkk�
%
% $Id: rts_smooth.m 109 2007-09-04 08:32:58Z jmjharti $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [M,S,D] = rts_smooth_sr(M,S,M_,S_,A)

  %
  % Check which arguments are there
  %
  if nargin < 5
    error('Too few arguments');
  end
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  S(:,:,end) = S(:,:,end)*S(:,:,end)'; 
  for k=(size(M,2)-1):-1:1
    
    D(:,:,k) = (S(:,:,k)*S(:,:,k)') * A(:,:,k)' / (S_(:,:,k)*S_(:,:,k)');
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - M_(:,k));
    S(:,:,k) = S(:,:,k)*S(:,:,k)' + D(:,:,k) * (S(:,:,k+1) - (S_(:,:,k)*S_(:,:,k)')) * D(:,:,k)';
  end

