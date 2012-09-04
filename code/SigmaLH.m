function [lh,glh] = SigmaLH(p,y,f,h,wi,ei,m0,P0)
    [~,~,~,~,~,lh,glh] = SigmaFilter(p,y,f,h,wi,ei,m0,P0);
    lh = -1*lh;
    glh = -1*glh;




