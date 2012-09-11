function [M,P,D] = rtssr_smooth(M,P,M_,P_,A)
	
  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end

  %
  % Extend A and Q if they are NxN matrices
  %
  if size(A,3)==1
    A = repmat(A,[1 1 size(M,2)]);
  end
%   if size(Q,3)==1
%     Q = repmat(Q,[1 1 size(M,2)]);
%   end

  %
  % Run the smoother
  %
  D = zeros(size(M,1),size(M,1),size(M,2));
  P(:,:,end) = P(:,:,end) * P(:,:,end)';
  for k=(size(M,2)-1):-1:1
    P(:,:,k) = P(:,:,k) * P(:,:,k)';
    P_(:,:,k) = P_(:,:,k) * P_(:,:,k)';
    
    %P_pred   = A(:,:,k) * P(:,:,k) * A(:,:,k)' + Q(:,:,k);
    D(:,:,k) = P(:,:,k)*A(:,:,k)' / P_(:,:,k); % smoother gain
    M(:,k)   = M(:,k) + D(:,:,k) * (M(:,k+1) - M_(:,k));
    P(:,:,k) = P(:,:,k) + D(:,:,k) * (P(:,:,k+1) - P_(:,:,k)) * D(:,:,k)';
    
%     Pj = [P(:,:,k),             D(:,:,k)*P(:,:,k+1);
%           P(:,:,k+1)*D(:,:,k)', P(:,:,k+1)];
      
    D(:,:,k) = D(:,:,k)*P(:,:,k+1); % the cross timestep covariance matrix
    %Pj =blkdiag(P(:,:,k), P(:,:,k+1));
%     try 
%         S = chol(P(:,:,k),'lower');
%         S = chol(P(:,:,k+1),'lower');
%         S = chol(Pj,'lower');
%     catch err
%         disp(k);
%         disp(Pj);
%         Ps = real(sqrtm(Pj));
%         Ps
%         Ps*Ps'
%         rethrow(err);
%      end
    
  end