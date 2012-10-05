function tmp=HarmonicRealInnerLH(QW,QX,Nqw)
%global dt c m0 P0 h f Jh Jf
load('../data/Harmonic_lqw_lqx_10_3751.mat');
tmp = zeros(Nqw,1);
p = p0;
for j=1:Nqw   
    p(gi) = [QW(j) QX(j)];
    tmp(j) = Harmonic_LH(p,Y);
end








