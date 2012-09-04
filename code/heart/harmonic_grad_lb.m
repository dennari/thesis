function [dSig,dmu,dQ,df,dR,dh] = harmonic_grad_lb(theta,w,dt,Qw)
%harmonic_grad_lb Gradient of the lower bound (Q-function) for the
% harmonic model
% theta = int, 1=Qx,2=Qw

%                         +-                                        -+
%                         |         sin(2 dt w)               2      |
%                         |  dt w - -----------,     sin(dt w)       |
%                         |              2                           |
% dQ/dQx =(1/(2*w)) *     |                                          |
%                         |               2      sin(2 dt w)         |
%                         |      sin(dt w) ,     ----------- + dt w  |
%                         |                           2              |
%                         +-                                        -+ 

%           +-         -+
%           |  0, 0, 0  |
%           |           |
% dQ/dQw =  |  0, 0, 0  |
%           |           |
%           |  0, 0, 1  |
%           +-         -+
  
  if nargin < 4
      Qw = 0;
  end
  zer = num2cell(zeros(1,6));
  [dSig,dmu,dQ,df,dR,dh] = deal(zer{:});
  
  if theta == 1 % dQ/dQx
    dQ = sinusoid_Q(w,1,Qw,dt); 
  elseif theta == 2 % dQ/dQw
    dQ = zeros(3);
    dQ(3,3) = 1;
  end
    
   

  
end