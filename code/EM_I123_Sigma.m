function [I1,I2,I3] = EM_I123_Sigma(f,h,m0,Y,MS,PS,DD)
	


N = size(Y,2);
xDim = size(MS,1);

[eij,wij] = CKFPoints(2*xDim);
M = size(eij,2);
w = diag(wij);

d = MS(:,1)-m0;
I1 = PS(:,:,1)+d*d';
I2 = 0;
I3 = 0;	

older = 1:xDim;
newer = (xDim+1):(2*xDim);


for k=2:N
    
    j = [k k+1];
    
    [xk1,xk] = jsig(MS(:,j),PS(:,:,j),DD(:,:,k));
    
	d = xk-f(xk1);
	I2 = I2 + d*w*d';
	
	d = repmat(Y(:,k),1,M)-h(xk);
    I3 = I3 + d*w*d';

end

function [old,new]=jsig(m,P,D)
    mj = m(:);
    Pj = [P(:,:,1) D; D'  P(:,:,2)];
    mr = repmat(mj,1,M); 
    try 
        %Ps = real(sqrtm(Pj));
        Ps = chol(Pj,'lower');
        xi = mr+Ps*eij;
    catch err
        disp(k);
        disp(Pj);
        rethrow(err);
    end
	% sigma points for x_k
    new = xi(newer,:);
    % sigma points for x_(k-1)
	old = xi(older,:);
end


end




