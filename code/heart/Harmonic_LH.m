function [lh,glh,varargout] = Harmonic_LH(p,y,gi,mult)
% parameters are 
% p(1)=qw,    angular velocity variance
% p(2)=r,     measurement variance
% p(3:3+c-1)  the signal component variances

  global m0 P0 H c
    
  if nargin < 4
    mult = 1;
  end
  if nargin < 3
    gi = [];
  end
    
  

  qx = p(3:end);
  if numel(p) < c+2
    qx = [qx repmat(p(end),1,c+2-numel(p))];
  end

  SQ = chol(sinusoid_Q(p(1),qx),'lower');
  SR = sqrt(p(2));
 
  h = @(x) H*x;
  f = @(x) sinusoid_f(x);

  xDim = size(m0,1);
  [usig,w] = CKFPoints(xDim);
  w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
  
  N = size(y,2)-1; % it is assumed that x0 has y0
  MM = zeros(xDim,N+1); MM(:,1) = m0;
  SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
  dm0 = zeros(xDim,1); dP0 = zeros(xDim);
  m = m0; S = P0; lh = 0; 
  glh = zeros(numel(gi),1); dm = dm0; dP = dP0;
  for k=1:(N+1)
      [m_,S_] = SigmaKF_Predict(m,S,f,SQ,[],usig,w);
      if k==N+1; break; end; 

      [m,S,~,IM,IS] = SigmaKF_Update(m_,S_,y(:,k+1),h,SR,usig,w);
      MM(:,k+1) = m;
      SS(:,:,k+1) = S;
      lh = lh + likelihood(y(:,k+1)-IM,IS*IS');

      if ~isempty(gi)
        [dm,dP,GLH] = Harmonic_Sens(dm,dP,A,H,gi,y(:,k+1)-IM,...
                        IS*IS',S_*S_'*H', MM(:,k),SS(:,:,k)*SS(:,:,k)');
        glh = glh + GLH;
      end
  end
  lh = mult*lh;
  glh = mult*glh;
  if nargout > 2
    varargout{1} = MM;
    varargout{2} = SS;
  end

end

