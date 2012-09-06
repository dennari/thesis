function [mF,SF] = SigmaSmootherSR(f,SQ,ms,Ss,m0,S0,wi,ei)
	
	XD = size(ms,1);
	N = size(ms,2);
	
    if nargin < 7 % use CKF by default
		[ei,wi] = CKFPoints(XD);
    end
	
    
	% the number of sigma points
	NS = size(ei,2);
	% it's possible to give different weights for the mean and variance
	if size(wi,1)==1
		wi = [wi;wi];
    end
    wm = wi(1,:)'; % weights for the mean
    wp = eye(NS)/sqrt(NS);%diag(wi(2,:)); % weights for the covariance
    
    
    
	%%%%%%%%%%%%%%%
	% SMOOTHER %%%%
	%%%%%%%%%%%%%%%
	SF = zeros(2*XD,2*XD,N);
	mF = zeros(2*XD,N);
	
	m_s = ms(:,end);
	S_s = Ss(:,:,end);

    
    for k=N-1:-1:0
	    if k > 0
	    	m = ms(:,k);
	    	S = Ss(:,:,k);
	    else
	    	m = m0;
	    	S = S0;
        end
        
        mr = repmat(m,1,NS); 
        sig = mr+S*ei;
        
        % propagate through the dynamics function
        sig_=f(sig);
        % apply the integration formula, wi(1,:) are the weights for mean
        m_ = sig_*wm;
        mr_ = repmat(m_,1,NS);

        
       
        X = (sig-mr)*wp;
        X_ = (sig_-mr_)*wp;
        
        [~,U] = qr([X_ SQ; X zeros(XD,XD)]',0);
        U = U';
        

	    % store the current values
	    m_s_ = m_s;
	    S_s_ = S_s;
	    

	    % smoother gain
	    G = U((XD+1):end,1:XD)/U(1:XD,1:XD);
	    
        m_s = m + G*(m_s-m_);
	    
        
        [~,S_s] = qr([U((XD+1):end,(XD+1):end) G*S_s_]',0);
        S_s = S_s';
        
	    % remember m_s_ and P_s_ have one higher k!
        % with this order m3(1:xDim,:) gets m_1:T - m_T:T
        % and m_0:T is m3(xDim+1:end,1)
        P_s_ = (S_s_*S_s_');
        P_s = S_s*S_s';
        %D = G*P_s_; % sqrt of the cross-covariance matrix
        %m3 = [m_s_;m_s];
	    %S3 = [Pp D; D' P];
        
        
        m3 = [m_s;m_s_];
        S3 = [P_s G*P_s_; P_s_*G' P_s_];
        
        %[I2_,I3_] = LBIntegrals(m3,P3);
        %I2 = I2+I2_;I3 = I3+I3_;
        
	    mF(:,k+1) = m3;
	    SF(:,:,k+1) = S3;
        
    end


    
end