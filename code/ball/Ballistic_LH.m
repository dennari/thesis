function [lh,glh,varargout] = Ballistic_LH(p,y,gi)
% parameters are 
%
% p{1}=lqx,  x process variance
% p{2}=lqy,  y process variance
% p{3}=lr,   measurement variance

    global A H P0 u m0
    if nargin < 4
      mult = 1;
    end
    if nargin < 3
      gi = [];
    end
    

    Q = ballisticQ2D(p(1),p(2));
    R = ballisticR(p(3:end));

    
    xDim = size(A,1);
    N = size(y,2)-1; % it is assumed that x0 has y0
    MM = zeros(xDim,N+1); MM(:,1) = m0;
    MM_ = zeros(xDim,N+1);
    PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
    PP_ = zeros(xDim,xDim,N+1);
    dm0 = zeros(xDim,1); dP0 = zeros(xDim);
    m = m0; P = P0; lh = 0; 
    glh = zeros(numel(gi),1); dm = dm0; dP = dP0;
    for k=1:(N+1)
        [m_,P_] = kf_predict(m,P,A,Q,[],u);
        MM_(:,k) = m_;
        PP_(:,:,k) = P_;
        if k==N+1; break; end; 
        
        [m,P,~,IM,IS] = kf_update(m_,P_,y(:,k+1),H,R);
        MM(:,k+1) = m;
        PP(:,:,k+1) = P;
        lh = lh + likelihood(y(:,k+1)-IM,IS);
        
        if ~isempty(gi)
          [dm,dP,GLH] = ballistic_sensitivityeq(dm,dP,A,H,gi,y(:,k+1)-IM,...
                          IS,P_*H', MM(:,k),PP(:,:,k),p);
          glh = glh + GLH;
        end
    end
    if nargout > 2
      varargout{1} = MM;
      varargout{2} = PP;
      varargout{3} = MM_;
      varargout{4} = PP_;
    end

end

