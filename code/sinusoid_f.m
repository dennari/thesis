function x_ = sinusoid_f(x,k,p,dt)
	A = @(w) [(cos(dt*w)), (sin(dt*w)/w),   0;
          (-w*sin(dt*w)), (cos(dt*w)),  0;
          0,                0,          1.001];
	x_ = zeros(size(x));
    for k=1:size(x,2)
       x_(:,k) = A(max(abs(x(3,k)),1e-9))*x(:,k); 
    end
	

end