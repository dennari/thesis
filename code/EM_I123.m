function [I1,I2,I3] = EM_I123(A,H,m0,Y,MS,PS,DD)
% DD - smoother gains

N = size(Y,2)-1; % it is assumed that there's a 'measurement' for x0
d = MS(:,1)-m0;
I1 = PS(:,:,1)+d*d';
I2 = 0;
I3 = 0;	

% I2o = 0;
% I3o = 0;	

xDim = size(A,1);
yDim = size(H,1);
%older = 1:xDim;
%newer = (xDim+1):(2*xDim);
%X00 = 0;
%X10 = 0;
%X11 = 0;

for k=1:N
    %j = [k k+1];
    %mj = MS(:,j);
    %Pj = PS(:,:,j);
    %D = DD(:,:,k);
    y = Y(:,k+1);
    
    mj = [MS(:,k);MS(:,k+1)];
    Pj = [PS(:,:,k)   DD(:,:,k)*PS(:,:,k+1); 
          PS(:,:,k+1)*DD(:,:,k)'  PS(:,:,k+1)];
%     
     I2 = I2 + Pj + mj*mj';
     X00 = PS(:,:,k+1)+MS(:,k+1)*MS(:,k+1)';
%     
    xyj = [X00          MS(:,k+1)*y';
           y*MS(:,k+1)' y*y'];
    I3 = I3 + xyj;
    
%     X11 = PS(:,:,k)+MS(:,k)*MS(:,k)';
%     X10 = DD(:,:,k)*PS(:,:,k+1)+MS(:,k)*MS(:,k+1)';
     
%     YY = y*y';
%     XY = MS(:,k+1)*y';
% 
%      I2o = I2o + X00 - A*X10 - (A*X10)' + A*X11*A';
%      I3o = I3o + YY - H*XY - (H*XY)' + H*X00*H';
        
end
%disp('----I123--------');

%glb = Q\(I2(older,newer)-A*I2(older,older));

IA = [A -eye(xDim)]';
IH = [H -eye(yDim)]';

I2 = IA'*I2*IA;
%disp(I2);
I3 = IH'*I3*IH;
%I2 = I2o;
%I3 = I3o;
%disp(I2-I2o);
%disp(I3);
%disp(I3o);

% I1 =   log(det(P0))+trace(I1/P0);
% I2 = N*log(det(Q))+trace(I2/Q);
% I3 = N*log(det(R))+trace(I3/R);
% 
% lb = -0.5*(I1+I2+I3);


end
