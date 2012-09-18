function [I1,I2,I3] = EM_I123(A,H,m0,Y,MS,PS,DD)
% DD - smoother gains

N = size(Y,2)-1; % it is assumed that there's a 'measurement' for x0
d = MS(:,1)-m0;
I1 = PS(:,:,1)+d*d';
I2 = 0;
I3 = 0;	

xDim = size(A,1);
yDim = size(H,1);

for k=1:N

  y = Y(:,k+1);

  mj = [MS(:,k);MS(:,k+1)];
  Pj = [PS(:,:,k)               DD(:,:,k)*PS(:,:,k+1); 
        PS(:,:,k+1)*DD(:,:,k)'  PS(:,:,k+1)];

  I2 = I2 + Pj + mj*mj';
  X00 = PS(:,:,k+1)+MS(:,k+1)*MS(:,k+1)';
  xyj = [X00          MS(:,k+1)*y';
         y*MS(:,k+1)' y*y'];
  I3 = I3 + xyj;
    
        
end


IA = [A -eye(xDim)]';
IH = [H -eye(yDim)]';

I2 = IA'*I2*IA;
I3 = IH'*I3*IH;



end
