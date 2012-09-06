function [ms,Ss,lh] = SigmaFilterSR(p,y,f,h,wi,ei,m0,S0)

	% f = model.f;
	% h = model.h;
	% wi = model.wi;
	% ei = model.ei;
	% m0 = model.m0;
	% P0 = model.P0;

	N = size(y,2);

	
	SQ = p{2}; % square root of Q
	SR = p{4}; % square root of R


	XD = size(p{1},1);
	YD = size(y,1);
    
    if isempty(wi) % use CKF by default
		[ei,wi] = CKFPoints(XD);
    end
	if isempty(f)
		% default to linear, A = p{1}
		f = @(x,k,p) p{1}*x;
	end
	if isempty(h)
		% default to linear, H = p{3}
		h= @(x,k,p) p{3}*x;
	end    
    
	% the number of sigma points
	NS = size(ei,2);
	% it's possible to give different weights for the mean and variance
	if size(wi,1)==1
		wi = [wi;wi];
    end
    wm = wi(1,:)'; % weights for the mean
    wp = eye(NS)/sqrt(NS);%diag(wi(2,:)); % weights for the covariance
    
	ms = zeros(XD,N);
	Ss = zeros(XD,XD,N);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FILTER  %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	m = m0;
	S = S0;
    lh = 0;
    

	for k=1:N
        % PREDICTION, p(x_k|y_1:k-1)=N(mean_pred,P_)
        
        % affine transform the unit sigma points
        % sig_,sigp_
        
        mr = repmat(m,1,NS); 
        sig = mr+S*ei;

        
        % propagate through the dynamics function
        sig_=f(sig,k,p);
        % apply the integration formula, wi(1,:) are the weights for mean
        m_ = sig_*wm;
        mr_ = repmat(m_,1,NS);
        
        % apply the integration formula, wi(2,:) are the weights for covariance
        [~,S_] = qr([(sig_-mr_)*wp SQ]',0);
        S_ = S_'; % we want lower triangular
        
%         if k < 5
%             P_ = S_*S_';
%             P_ = P_(:);
%             P = S*S';
%             P = P(:);
%             
%             fprintf(1,'sqrt ( %.2f, %.2f, %.2f, %.2f )\n',m'*m,P'*P,m_'*m_,P_'*P_); 
%         end
        
        %ms_(:,k) = m_;
        %Ss_(:,:,k) = S_;
        
        % affine transform the unit sigma points
        sigu = mr_+S_*ei;
        % propagate through the measurement function
        sigu_ = h(sigu,k,p);
        % integrate
        y_ = sigu_*wm;
        yr_ = repmat(y_,1,NS);
        
        % create a big matrix
        X = (sigu-mr_)*wp;
        Z = (sigu_-yr_)*wp;
        [~,T] = qr([Z SR; X zeros(XD,YD)]',0);
        T = T';
    
        Sy = T(1:YD,1:YD);
        C = T((YD+1):end,1:YD);
        
        % the gain
        K = C/Sy;
        
        
        d = y(:,k)-y_;
        m = m_+K*d;
        S = T((YD+1):end,(YD+1):end);
        ms(:,k) = m;
        Ss(:,:,k) = S;
        
        Py = Sy*Sy'; % innovations covariance matrix
		lh = lh - 0.5*(d'/Py*d + log(det(Py)) + YD*log(2*pi));
        
        %mloglik - log(2.*pi).*(no/2)- log(det(Sy*Sy'))/2 - resid'/(Sy*Sy')*resid/2;
        

	end

end




