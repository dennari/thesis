function [lh,glh,varargout] = Ballistic_LH(p,y,gi)
% parameters are 
%
% p{1}=lqx,  x process variance
% p{2}=lqy,  y process variance
% p{3}=lr,   measurement variance
% p{5}=ux,    constant input x-component
% p{6}=uy,    constant input y-component 

    global A H P0 m0

    if nargin < 3
      gi = [];
    end
    

    Q = ballisticQ2D(p(1),p(2));
    R = ballisticR(p(3));
    u = ballisticU(p(5),p(6));

    
    xDim = size(A,1);
    N = size(y,2)-1; % it is assumed that x0 has y0
    MM = zeros(xDim,N+1); MM(:,1) = m0;
    MM_ = zeros(xDim,N+1);
    PP = zeros(xDim,xDim,N+1); PP(:,:,1) = P0;
    PP_ = zeros(xDim,xDim,N+1);
    dm0 = zeros(xDim,1); dP0 = zeros(xDim);
    m = m0; P = P0; lh = 0; 
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
        [m_,P_] = kf_predict(m,P,A,Q,[],u);
        MM_(:,k) = m_;
        PP_(:,:,k) = P_;
        
        % run the partial derivative predictions
        for i=1:numel(gi)
          [dm,dP] = dd{i}{:};
          [dm_,dP_] = ballistic_dpred(gi(i),dm,dP,m,P,p);
          %[dm_,dP_] = dSigmaKF_Predict(m,m_,S,f,dm,dP,dQ,Jf,usig,w);
          dd_{i} = {dm_,dP_}; 
        
          
        end
        
        if k==N+1; break; end; 
        
        [m,P,K,my,Sy] = kf_update(m_,P_,y(:,k+1),H,R);
        MM(:,k+1) = m;
        PP(:,:,k+1) = P;
        lh = lh + likelihood(y(:,k+1)-my,Sy);
        
        % run the partial derivative updates
        for i=1:numel(gi)
          [dm_,dP_] = dd_{i}{:};
          %[dm,dP,dmy,dSy] = dSigmaKF_Update(m_,S_,h,dm_,dP_,K,my,Sy,yy,dR,Jh,usig,w);
          [dm,dP,dmy,dSy] = ballistic_dupdate(gi(i),dm_,dP_,P_,K,my,Sy,y(:,k+1),p);
          dd{i} = {dm,dP};
          glh(i) = glh(i) + dlikelihood(Sy,dSy,y(:,k+1),my,dmy);
        end

    end
    if nargout > 2
      varargout{1} = MM;
      varargout{2} = PP;
      varargout{3} = MM_;
      varargout{4} = PP_;
      varargout{5} = Q;
      varargout{6} = u;
    end

end

