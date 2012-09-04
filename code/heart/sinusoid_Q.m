function Q_ = sinusoid_Q(w,Qc,Qw,dt)
% The signal part of Q is    
%           +-                                        -+
%           |                                   2      |
%           |  dt w   sin(2 dt w)      sin(dt w)       |
%           |  ---- - -----------,     ----------      |
%           |   2          4               2           |
% (Qx/w) *  |                                          |
%           |               2                          |
%           |      sin(dt w)       sin(2 dt w)   dt w  |
%           |      ----------,     ----------- + ----  |
%           |          2                4         2    |
%           +-                                        -+    
        
    Q11 = (2*dt*w - sin(2*dt*w))/2;
    Q12 = sin(dt*w)^2;
    Q22 = (2*dt*w + sin(2*dt*w))/2;
    
    Q_ = blkdiag((Qc/(2*w))*[Q11 Q12; Q12 Q22],Qw);
   
end

%     Q11 = dt^2/3;
%     Q12 = dt/2;
%     Q22 = 1;