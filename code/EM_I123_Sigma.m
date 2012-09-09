function [I1,I2,I3] = EM_I123_Sigma(f,h,m0,Y,MS,PS,DD)
	


N = size(Y,2);
xDim = size(MS,1)/2;
yDim = size(Y,1);

[eij,wij] = CKFPoints(2*xDim);
M = size(eij,2);
% weight matrix
wmx = repmat(wij(1,:),xDim,1);
wmy = repmat(wij(1,:),yDim,1);

d = MS(:,1)-m0;
I1 = PS(:,:,1)+d*d';
I2 = 0;
I3 = 0;	

newer = 1:xDim;
older = (xDim+1):(2*xDim);


for k=1:N
	% compute integrals I_2 and I_3 that are w.r.t N(m3,P3)
    
    j = [k k+1];
    % get the sigma points from the joint distribution
    [xk1,xk] = jsig(MS(:,j),PS(:,:,j),DD(:,:,k));
    
	d = xk-f(xk1);
	I2 = I2 + (wmx.*d)*d';
	
	d = repmat(Y(:,k),1,M)-h(xk);
    I3 = I3 + (wmy.*d)*d';

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




