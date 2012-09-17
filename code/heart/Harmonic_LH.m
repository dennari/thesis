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
  Jh = @(x) H;
  f = @(x) sinusoid_f(x);
  Jf = @(x) sinusoid_Jf(x);
  
  xDim = size(m0,1);
  [usig,w] = CKFPoints(xDim);
  w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
  
  N = size(y,2)-1; % it is assumed that x0 has y0
  MM = zeros(xDim,N+1); MM(:,1) = m0;
  SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
  dm0 = zeros(xDim,1); dP0 = zeros(xDim);
  m = m0; S = P0; lh = 0; 
  glh = zeros(numel(gi),1);
  
  
  % initialize the partial derivative structures
  if ~isempty(gi)
    dd = cell(1,numel(gi)); 
    for i=1:numel(gi)
      dd{i} = {dm0,dP0};
    end
    dd_ = cell(1,numel(gi));
  end
  
  
  for k=1:(N+1)
      [m_,S_] = SigmaKF_Predict(m,S,f,SQ,[],usig,w);
      
      % run the partial derivative predictions
      for i=1:numel(gi)
        [dm,dP] = dd{i}{:};
        [dQ,~] = dQdR(gi(i));
        [dm_,dP_] = dSigmaKF_Predict(m,m_,S,f,dm,dP,dQ,Jf,usig,w);
        dd_{i} = {dm_,dP_}; 
      end
      
      if k==N+1; break; end; 
      
      yy = y(:,k+1);
      [m,S,K,my,CSy,CC] = SigmaKF_Update(m_,S_,yy,h,SR,usig,w);
      %%% CSy and CC are Cholesky decompositions _HERE_ %%%
      Sy = CSy*CSy';
      C = CC*CC';
      MM(:,k+1) = m;
      SS(:,:,k+1) = S;
      lh = lh + likelihood(yy-my,Sy);
      
      % run the partial derivative updates
      for i=1:numel(gi)
        [dm_,dP_] = dd_{i}{:};
        [~,dR] = dQdR(gi(i));
        [dm,dP,dmy,dSy] = dSigmaKF_Update(m_,S_,h,dm_,dP_,K,my,Sy,yy,dR,Jh,usig,w);
        dd{i} = {dm,dP};
        glh(i) = glh(i) + dlikelihood(Sy,dSy,yy,my,dmy);
      end
      
  end
  lh = mult*lh;
  glh = mult*glh;
  if nargout > 2
    varargout{1} = MM;
    varargout{2} = SS;
  end
  
end

function [dQ,dR]=dQdR(i)
  global P0 H c
  
  dQ = zeros(size(P0));
  dR = zeros(size(H,1));

  if(i==1) % dlb/dqw
      dQ(1,1) = 1;
  end
  if(i >= 3) % dlb/dqx(ri)
      %dQ = sinusoid_Q(0,dqxi(ri,:),dt);
      %ri = ri + 1;
      %if sum(gi>=3) > 1
      %  wh = zeros(1,c);
      %  wh(gi-2) = 1;
      %  dQ = sinusoid_Q(0,wh);
      %else  
        dQ = sinusoid_Q(0,ones(1,c));
      %end
  end
  if(i==2) % dlb/dr
      dR = 1;
  end
  
end

