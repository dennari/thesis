function [lh,glh,varargout] = Ballistic_LH_SR(p,y,gi)
% parameters are 
%
% p{1}=lqx,  x process variance
% p{2}=lqy,  y process variance
% p{3}=lr,   measurement variance

    global A H P0 u m0

    if nargin < 3
      gi = [];
    end
    

    Q = ballisticQ2D(p(1),p(2));
    R = ballisticR(p(3:end));

    SQ = chol(Q,'lower');
    SR = chol(R,'lower');
    
    xDim = size(A,1);
    N = size(y,2)-1; % it is assumed that x0 has y0
    MM = zeros(xDim,N+1); MM(:,1) = m0;
    MM_ = zeros(xDim,N+1);
    SS = zeros(xDim,xDim,N+1); SS(:,:,1) = P0;
    SS_ = zeros(xDim,xDim,N+1);
    dm0 = zeros(xDim,1); dP0 = zeros(xDim);
    m = m0; S = P0; lh = 0; 
    glh = zeros(numel(gi),1); dm = dm0; dP = dP0;
    for k=1:(N+1)
        [m_,S_] = kfsr_predict(m,S,A,SQ,[],u);
        MM_(:,k) = m_;
        SS_(:,:,k) = S_;
        if k==N+1; break; end; 
        
        [m,S,~,my,CSy] = kfsr_update(m_,S_,y(:,k+1),H,SR);
        Sy = CSy*CSy';
        MM(:,k+1) = m;
        SS(:,:,k+1) = S;
        lh = lh + likelihood(y(:,k+1)-my,Sy);
        
        if ~isempty(gi)
          [dm,dP,GLH] = ballistic_sensitivityeq(dm,dP,A,H,gi,y(:,k+1)-my,...
                          Sy,S_*S_'*H', MM(:,k),SS(:,:,k)*SS(:,:,k)',p);
          glh = glh + GLH;
        end
    end
    if nargout > 2
      varargout{1} = MM;
      varargout{2} = SS;
      varargout{3} = MM_;
      varargout{4} = SS_;
    end

end

