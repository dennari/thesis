function tmp=HarmonicRealInnerLB(QW,QX,Nqw,N,p0,P0,I1,I2,I3)

tmp = zeros(Nqw,1);

R = sinusoid_R(p0(2));
p = p0;
gi = [1 3];
for j=1:Nqw   
    p(gi) = [QW(j) QX(j)];
    Q = sinusoid_Q(p(1),p(3));
    Q(2,2) = 1; Q(4,4) = 1; Q(6,6) = 1;
    II1 =   log(det(P0))+trace(I1/P0);
    II2 = N*log(det(Q))+trace(I2/Q);
    II3 = N*log(det(R))+trace(I3/R);
    tmp(j) = -0.5*(II1+II2+II3);
    
end








