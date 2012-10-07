function [I1,I2,I3] = EM_I123_Sigma(f,h,m0,Y,MS,SS,DD)
	


N = size(Y,2)-1;  % it is assumed that there's a 'measurement' for x0
xDim = size(MS,1);

[eij,wij] = CKFPoints(2*xDim);
M = size(eij,2);
w = diag(wij);

d = MS(:,1)-m0;
I1 = SS(:,:,1)*SS(:,:,1)'+d*d';
I2 = 0;
I3 = 0;	

older = 1:xDim;
newer = (xDim+1):(2*xDim);


for k=1:N
  y = Y(:,k+1);
  j = [k k+1];

  [xk1,xk] = jsig(MS(:,j),SS(:,:,j),DD(:,:,k));

  d = xk-f(xk1);
  I2 = I2 + d*w*d';

  d = repmat(y,1,M)-h(xk);
  I3 = I3 + d*w*d';

end

function [old,new]=jsig(m,S,D)
    P1 = S(:,:,1)*S(:,:,1)';
    %P1(abs(P1) < 1e-8) = 0;
    P2 = S(:,:,2)*S(:,:,2)';
    %P2(abs(P2) < 1e-8) = 0;
    %D(abs(D) < 1e-8) = 0;
    
    %Pj = [P1 D*P2; P2*D'  P2];
    try 
        %Ps = chol(Pj,'lower');
        % Schur complement
        SC = P2-P2*D'/P1*D*P2;
        % force positive definiteness
        SC(2,2) = 1; SC(4,4) = 1; SC(6,6) = 1;
        CSC = chol(SC,'lower');
        CSC(2,2) = 0; CSC(4,4) = 0; CSC(6,6) = 0;
        
        S = [S(:,:,1)         zeros(size(S(:,:,1)));
             P2*D'/S(:,:,1)'  CSC];
        %S(abs(S) < 1e-12) = 0;
        %sum(S(:)-Ps(:))
        
        xi = repmat(m(:),1,M)+S*eij;
    catch err
        disp(k);
        %SC = P2-P2*D'/P1*D*P2;
        %disp(P1);
        %disp(P2);
        %disp(D);
        rethrow(err);
    end
    % sigma points for x_k
    new = xi(newer,:);
    % sigma points for x_(k-1)
    old = xi(older,:);
end


end




