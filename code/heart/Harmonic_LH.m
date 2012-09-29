function [lh,glh,varargout] = Harmonic_LH(p,y,gi)
% parameters are 
% p(1)=lqw,    log(sqrt) angular velocity variance
% p(2)=lr,     log(sqrt) measurement variance
% p(3:3+c-1)   log(sqrt) component variances

  global m0 P0 f h Jf Jh
    

  if nargin < 3
    gi = [];
  end
  
  Q = sinusoid_Q(p(1),p(3));
  if Q == diag(diag(Q))
    SQ = sqrt(Q);
  else
    SQ = chol(Q,'lower');
  end
  
  SR = sqrt(sinusoid_R(p(2)));
 
  
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
      [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);
      
      % run the partial derivative predictions
      for i=1:numel(gi)
        [dm,dP] = dd{i}{:};
        [dQ,~] = dQdR(gi(i),p);
        [dm_,dP_] = dSigmaKF_Predict(m,m_,S,f,dm,dP,dQ,Jf,usig,w);
        dd_{i} = {dm_,dP_}; 
      end
      
      if k==N+1; break; end; 
      
      yy = y(:,k+1);
      [m,S,K,my,CSy] = SigmaKF_Update(m_,S_,yy,h,SR,usig,w);
      %%% CSy and CC are Cholesky decompositions _HERE_ %%%
      Sy = CSy*CSy';
      MM(:,k+1) = m;
      SS(:,:,k+1) = S;
      lh = lh + likelihood(yy-my,Sy);
      
      % run the partial derivative updates
      for i=1:numel(gi)
        [dm_,dP_] = dd_{i}{:};
        [~,dR] = dQdR(gi(i),p);
        [dm,dP,dmy,dSy] = dSigmaKF_Update(m_,S_,h,dm_,dP_,K,my,Sy,yy,dR,Jh,usig,w);
        dd{i} = {dm,dP};
        glh(i) = glh(i) + dlikelihood(Sy,dSy,yy,my,dmy);
      end
      
  end
  if nargout > 2
    varargout{1} = MM;
    varargout{2} = SS;
  end
  if nargout > 4
    varargout{3} = SQ;
  end
  
end

function [dQ,dR]=dQdR(i,p)
  global c
  
  dQ = zeros(2*c+1);
  dR = 0;

  if(i==1) % dlh/dlog(qw)
      dQ(1,1) = 1*2*exp(2*p(1));
  end
  if(i >= 3) % dlb/dlog(qx(ri))
      %dQ = sinusoid_Q(0,dqxi(ri,:),dt);
      %ri = ri + 1;
      %if sum(gi>=3) > 1
      %  wh = zeros(1,c);
      %  wh(gi-2) = 1;
      %  dQ = sinusoid_Q(0,wh);
      %else  
        dQ = sinusoid_Q(0,1,1)*2*exp(2*p(3));
      %end
  end
  if(i==2) % dlb/dlog(r)
      dR = 1*2*exp(2*p(2));
  end
  
end

