function [mF,PF] = SigmaSmoother(ms,Ps,ms_,Ps_,Ds,m0,P0)
	
	XD = size(ms,1);
	N = size(ms,2);
	
	%%%%%%%%%%%%%%%
	% SMOOTHER %%%%
	%%%%%%%%%%%%%%%
	PF = zeros(2*XD,2*XD,N);
	mF = zeros(2*XD,N);
	
	m_s = ms(:,end);
	P_s = Ps(:,:,end);

    
    for k=N-1:-1:0
	    if k > 0
	    	m = ms(:,k);
	    	P = Ps(:,:,k);
	    else
	    	m = m0;
	    	P = P0;
	    end

	    % prediction
	    m_ = ms_(:,k+1);
	    P_ = Ps_(:,:,k+1);
	    
	    % the cross-covariance
	    C = Ds(:,:,k+1);
	    

	    % store the current values
	    m_s_ = m_s;
	    P_s_ = P_s;
	    

	    % update
	    G = C / P_;
	    m_s = m + G*(m_s-m_);
	    P_s = P +G*(P_s-P_)*G';
	    D = G*P_s_;
        
	    % remember m_s_ and P_s_ have one higher k!
        % with this order m3(1:xDim,:) gets m_1:T - m_T:T
        % and m_0:T is m3(xDim+1:end,1)
        %m3 = [m_s_;m_s];
	    %P3 = [P_s_ D; D P_s];
        
        
        m3 = [m_s;m_s_];
        P3 = [P_s G*P_s_; P_s_*G' P_s_];
        
        %[I2_,I3_] = LBIntegrals(m3,P3);
        %I2 = I2+I2_;I3 = I3+I3_;
        
	    mF(:,k+1) = m3;
	    PF(:,:,k+1) = P3;
        
    end


    
end