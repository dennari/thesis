function [glb] = EM_LB_AR2(vi,Q,R,P0,m0,m0T,N,I1,I2,I3)



% gradient

glb = zeros(numel(vi));
for j=1:numel(vi)
    % assume all zero
    dSig = zeros(size(P0)); dQ = zeros(size(Q));dR = zeros(size(R));
    dmu = zeros(size(Q,1),1);
    dA = dQ;

    if(vi(1)==1) % dlb/a1
        dA(1,1) = 1;
    end
    if(vi(j)==3) % dlb/dq1
        dQ(1,1) = 1;
    end
     if(vi(j)==5) % dlb/dr
        dR(1,1) = 1;
    end
   
    % x_0
    glb1 = P0\(dSig/P0*I1+2*dmu*(m0T-m0)'-dSig);
    % x_k|x_(k-1)
    glb2 = Q\(dQ/Q*I2-N*dQ);
    % y_k|x_k
    glb3 = R\(dR/R*I3-N*dR);

    glb(j) = 0.5*(trace(glb1)+trace(glb2)+trace(glb3));
    %glb(j) = 0.5*trace(glb2);
    
end

end




