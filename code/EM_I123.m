function [I1,I2,I3] = EM_I123(A,H,m0,Y,MS,PS,DD,u)
% DD - smoother gains

N = size(Y,2)-1; % it is assumed that there's a 'measurement' for x0
d = MS(:,1)-m0;
I1 = PS(:,:,1)+d*d';
I2 = 0;
I3 = 0;	




xDim = size(A,1);
yDim = size(H,1);

if nargin < 8
  u = zeros(xDim,1); % constant input
end
lr = xDim+1:2*xDim; 

for k=1:N

  y = Y(:,k+1);

  mj = [MS(:,k);MS(:,k+1)];
  Pj = [PS(:,:,k)               DD(:,:,k)*PS(:,:,k+1); 
        PS(:,:,k+1)*DD(:,:,k)'  PS(:,:,k+1)];

  XX = [Pj+mj*mj' mj*u'
        u*mj'     u*u'];
  
  I2 = I2 + XX;  
  
  YY = [XX(lr,lr)    MS(:,k+1)*y';
         y*MS(:,k+1)' y*y'];
  
  I3 = I3 + YY;
    
        
end


IA = [A -eye(xDim) eye(xDim)]';
IH = [H -eye(yDim)]';

I2 = IA'*I2*IA;
I3 = IH'*I3*IH;


end
