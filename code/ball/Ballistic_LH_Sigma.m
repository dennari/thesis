function [lh,glh,varargout] = Ballistic_LH_Sigma(p,y,gi,mult)
% parameters are 
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement _STD_

    global A H P0 u
    if nargin < 4
      mult = 1;
    end
    if nargin < 3
      gi = [];
    end
    

 
    SQ = chol(ballisticQ2D(p(3),p(4)),'lower');
    SR = sqrt(p(5))*eye(size(y,1));
    %m0 = [0 p(1) g0x 0 p(2) g0y]';
    m0 = [0 p(1) 0 p(2)]';
    f = @(x) A*x;
    h = @(x) H*x;
    Jf = @(x) A;
    Jh = @(x) H;
    
    xDim = size(A,1);
    [usig,w] = CKFPoints(size(A,1));
    w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
    N = size(y,2)-1; % it is assumed that x0 has y0
    MM = zeros(xDim,N+1); MM(:,1) = m0;
    SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
    dm0 = zeros(xDim,1); dP0 = zeros(xDim);
    m = m0; S = P0; lh = 0; 
    glh = zeros(numel(gi),1); dm_o = dm0; dP_o = dP0;
    
      % initialize the partial derivative structures
    if ~isempty(gi)
      dd = cell(1,numel(gi)); 
      for i=1:numel(gi)
        dd{i} = {dm0,dP0};
      end
      dd_ = cell(1,numel(gi));
    end
    
    
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,u,usig,w);
        
        % run the partial derivative predictions
        for i=1:numel(gi)
          [dm,dP] = dd{i}{:};
          [dQ,~] = dQdR(gi(i));
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
  global A H
  
    dA = zeros(size(A));
    dQ = zeros(size(A));
    dR = zeros(size(H,1));

    if(i==1) % dlh/dv0x
        dA(1,1) = 1;
    end
    if(i==2) % dlh/dv0y
        dA(1,2) = 1;
    end
    if(i==3) % dlh/dqx
        dQ = ballisticQ2D(1,0);
    end
    if(i==4) % dlh/dqy
        dQ = ballisticQ2D(0,1);
    end
    if(i==5) % dlh/dr
        dR = eye(size(dR,1));
    end
  
end



