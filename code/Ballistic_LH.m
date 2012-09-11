function [lh,glh,varargout] = Ballistic_LH(p,y,gi,mult)
% parameters are 
% p{1}=v0x, x component of the mean of the initial velocity
% p{2}=v0y, y component of the mean of the initial velocity 
% p{3}=qx,  x process variance
% p{4}=qy,  y process variance
% p{5}=r,   measurement _STD_

    global A H P0 g0x g0y
    if nargin < 4
      mult = 1;
    end
    if nargin < 3
      gi = [];
    end
    
    if iscell(p)
      Q = p{1};
      R = p{2};
      m0 = [0 p{3} g0x 0 p{4} g0y]';
    else
      Q = ballisticQ(p(3),p(4));
      R = p(5)^2*eye(size(y,1));
      m0 = [0 p(1) g0x 0 p(2) g0y]';
    end
    
    
    xDim = size(A,1);
    N = size(y,2);
    MM = zeros(xDim,N+1); MM(:,1) = m0;
    PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
    dm0 = zeros(xDim,1); dP0 = zeros(xDim);
    m = m0; P = P0; lh = 0; 
    glh = zeros(numel(gi),1); dm = dm0; dP = dP0;
    for k=1:(N+1)
        [m_,P_] = kf_predict(m,P,A,Q);
        
        if k==N+1; break; end; 
        
        [m,P,~,IM,IS,LH] = kf_update(m_,P_,y(:,k),H,R);
        MM(:,k+1) = m;
        PP(:,:,k+1) = P;
        lh = lh + log(LH);
        
        if ~isempty(gi)
          [dm,dP,GLH] = ballistic_sensitivityeq(dm,dP,A,H,gi,y(:,k)-IM,...
                          IS,P_*H', MM(:,k),PP(:,:,k));
          glh = glh + GLH;
        end
    end
    lh = mult*lh;
    glh = mult*glh;
    if nargout > 2
      varargout{1} = MM;
      varargout{2} = PP;
    end

end

