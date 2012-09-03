function Q_ = sinusoid_Q(m,k,p,Qc,Qw,dt)
    
    Q =@(w) [-(Qc*(sin(2*dt*w) - 2*dt*w))/(4*w^3),...          
              (Qc*sin(dt*w)^2)/(2*w^2);
              (Qc*sin(dt*w)^2)/(2*w^2),... 
              (Qc*(sin(2*dt*w) + 2*dt*w))/(4*w)];
	
    Q_ = zeros(3);      
    Q_(1:2,1:2) = Q(max(abs(m(3)),1e-9));
    Q_(end,end) = Qw;

end