function [lb,glb] = EM_LB_Linear(model,p,Y,m3,P3)
	
A = p{1};
Q = p{2};
H = p{3};
R = p{4};

xDim = size(A,1);
yDim = size(H,1);
N = size(Y,2);

I1 = 0;
I2 = 0;
I3 = 0;	


m0 = model.m0;
P0 = model.P0;

newer = 1:xDim;
older = (xDim+1):(2*xDim);


for k=1:N
    %m = rot90(m3(:,k),2);
	%P = rot90(P3(:,:,k),2);
    m = m3(:,k);
	P = P3(:,:,k);

    
    if k == 1
        % p(x_0|y_1:N)
        d = m(older)-m0;
        I1 = P(older,older)+d*d';
    end
    
    I2 = I2 + P+m*m';
    
    y = Y(:,k);
    I3 = I3 + [y*y'         y*m(newer)';
               m(newer)*y'  m(newer)*m(newer)'];
    
    %X00 = X00 + P(older,older)+m(older)*m(older)';
    %X01 = X01 + P(older,newer)+m(older)*m(older)';
    %X00 = X00 + P(older,older)+m(older)*m(older)';

        
end

glb = Q\(I2(older,newer)-A*I2(older,older));

IA = [eye(xDim) -A]';
IH = [eye(yDim) -H]';
I2 = IA'*I2*IA;
I3 = IH'*I3*IH;

I1 =   log(det(P0))+trace(I1/P0);
I2 = N*log(det(Q))+trace(I2/Q);
I3 = N*log(det(R))+trace(I3/R);

lb = -0.5*(I1+I2+I3);


end
