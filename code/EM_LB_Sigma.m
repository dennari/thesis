function [lb,varargout] = EM_LB_Sigma(model,p,Y,JM,JP)
	


N = size(Y,2);
xDim = size(JM,1)/2;
yDim = size(Y,1);

[eij,wij] = CKFPoints(2*xDim);
M = size(eij,2);
% weight matrix
wmx = repmat(wij(1,:),xDim,1);
wmy = repmat(wij(1,:),yDim,1);

I1 = 0;
I2 = 0;
I3 = 0;	

newer = 1:xDim;
older = (xDim+1):(2*xDim);

m0 = model.m0;
P0 = model.P0;

%IIx = [eye(xDim);-eye(xDim)];
%IIy = [eye(yDim);-eye(yDim)];

for k=1:N
	% compute integrals I_2 and I_3 that are w.r.t N(m3,P3)
	m = JM(:,k);
	P = JP(:,:,k);
    
    if k == 1
        % p(x_0|y_1:N)
        d = m(older)-m0;
        I1 = P(older,older)+d*d';
    end
    
    % get the sigma points from the joint distribution
    [xk1,xk] = jsig(m,P);
    
	d = xk-model.f(xk1,k,p);
	I2 = I2 + (wmx.*d)*d';
	
    %j = [xk ; model.f(xk1,k,p)];
	%I2 = I2 + IIx'*(wmx.*j)*j'*IIx;
	
	d = repmat(Y(:,k),1,M)-model.h(xk,k,p);
    I3 = I3 + (wmy.*d)*d';

end


Q = p{2};
R = p{4};




if isfield(model,'gradDim') && model.gradDim > 0
	% gradient
    glb = zeros(model.gradDim,1);
    for j=1:model.gradDim
        [dSig,dmu,dQ,df,dR,dh] = model.grad(j,p,JM,JP);
        % x_0
        glb1 = P0\(dSig/P0*I1+2*dmu*(JM(older,1)-m0)-dSig);
        % x_k|x_(k-1)
        glb2 = Q\(dQ/Q*I2+2*df-N*dQ);
        %glb2 = I2/Q^2-N/Q;
        
        % y_k|x_k
        glb3 = R\(dR/R*I3+2*dh-N*dR);
		
        glb(j) = 0.5*(trace(glb1)+trace(glb2)+trace(glb3));
    end
	varargout{1} = glb;
end

I1 =  log(det(P0))+trace(I1/P0);
I2 = N*log(det(Q))+trace(I2/Q);
I3 = N*log(det(R))+trace(I3/R);

lb = -0.5*(I1+I2+I3);


function [old,new]=jsig(m,P)
    mr = repmat(m,1,M); 
	xi = mr+chol(P,'lower')*eij;
	% sigma points for x_k
    new = xi(newer,:);
    % sigma points for x_(k-1)
	old = xi(older,:);
end

end




