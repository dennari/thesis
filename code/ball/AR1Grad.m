function [dSig,dmu,dQ,df,dR,dh] = AR1Grad(theta,p,JM,JP)
%AR1GRAD Gradient for the AR1 SSM model in EM
% theta = int, 1=A,2=Q,3=H,4=R
% xk1 = [xDim x M], sigma points for x_(k-1)
% xk = [xDim x M], sigma points for x_k
% wj = [1 x M], weights
% p = struct, parameters of the model    
    xDim = size(JM,1)/2;
    newer = 1:xDim;
    older = (xDim+1):(2*xDim);    
    [eij,wij] = CKFPoints(2*xDim);
    M = size(eij,2);
    
    zer = num2cell(zeros(1,6));
    [dSig,dmu,dQ,df,dR,dh] = deal(zer{:});

    A = p{1};
    switch theta
        case 1 % A
            for k=1:size(JM,2)
                [xk1,xk] = jsig(JM(:,k),JP(:,:,k));
                df = df + xk1*diag(wij)*xk'-xk1*diag(wij)*xk1'*A';
            end
        
        case 2 % Q
            dQ = 1;
    end
    
    function [old,new]=jsig(m,P)
        mr = repmat(m,1,M); 
        xi = mr+chol(P,'lower')*eij;
        % sigma points for x_k
        new = xi(newer,:);
        % sigma points for x_(k-1)
        old = xi(older,:);
    end
end