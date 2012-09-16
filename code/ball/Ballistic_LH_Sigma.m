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
    
    xDim = size(A,1);
    [usig,w] = CKFPoints(size(A,1));
    w(:,2) = sqrt(w(:,1)); % add weights for square root filt/smooth
    N = size(y,2)-1; % it is assumed that x0 has y0
    MM = zeros(xDim,N+1); MM(:,1) = m0;
    SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
    dm0 = zeros(xDim,1); dP0 = zeros(xDim);
    m = m0; S = P0; lh = 0; 
    glh = zeros(numel(gi),1); dm = dm0; dP = dP0;
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,u,usig,w);
        if k==N+1; break; end; 
        
        [m,S,~,IM,IS] = SigmaKF_Update(m_,S_,y(:,k+1),h,SR,usig,w);
        MM(:,k+1) = m;
        SS(:,:,k+1) = S;
        lh = lh + likelihood(y(:,k+1)-IM,IS*IS');
        
        if ~isempty(gi)
          [dm,dP,GLH] = ballistic_sensitivityeq(dm,dP,A,H,gi,y(:,k+1)-IM,...
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

