function tmp=HarmonicRealInnerLH(QW,QX,Nqw,H,Y,SR,usig,w,m0,P0,dt)
%'fghfghf'
%pwd()
Qf = @(qw,qx) diag([qw^2, 0, qx^2, 0, qx^2, 0, qx^2]);
   f = @(x) sinusoid_f(x,dt);
   h = @(x) H*x;
tmp = zeros(Nqw,1);
  N = size(Y,2)-1;
  for j=1:Nqw   
    SQ = sqrt(Qf(QW(j),QX(j)));
    lh = 0; m = m0; S = P0;
    
    for k=1:(N+1)
        [m_,S_] = SigmaKF_Predict(m,S,f,SQ,usig,w);

        if k==N+1; break; end; 

      
        [m,S,~,my,CSy] = SigmaKF_Update(m_,S_,Y(:,k+1),h,SR,usig,w);
        %%% CSy and CC are Cholesky decompositions _HERE_ %%%
        Sy = CSy*CSy';
        lh = lh + likelihood(Y(:,k+1)-my,Sy);


    end
    j
    tmp(j) = lh;
    %fprintf('%.2f\n',100*((i-1)*Nqw+j)/(Nqx*Nqw));
  end








