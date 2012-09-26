function [ u ] = ballisticU( gx,gy)
global dt
    
  u = [0 dt*gx 0 dt*gy]';

end

