function [lb] = EM_LB(P0,Q,R,N,I1,I2,I3)
	

I1 =  log(det(P0))+trace(P0\I1);
I2 = N*log(det(Q))+trace(Q\I2);
I3 = N*log(det(R))+trace(R\I3);

%I1 =   log(det(P0))+trace(I1/P0);
%I2 = N*log(det(Q))+trace(I2/Q);
%I3 = N*log(det(R))+trace(I3/R);

lb = -0.5*(I1+I2+I3);


end




