function [ms,Ps,ms_,Ps_,Ds,varargout] = SigmaFilter(p,y,f,h,wi,ei,m0,P0)

	% f = model.f;
	% h = model.h;
	% wi = model.wi;
	% ei = model.ei;
	% m0 = model.m0;
	% P0 = model.P0;

	N = size(y,2);

	no = nargout - 5;
	
	Q = p{2};
	R = p{4};


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

	ms = zeros(XD,N);
	ms_ = zeros(XD,N);
	Ps = zeros(XD,XD,N);
	Ps_ = zeros(XD,XD,N);
	Ds = zeros(XD,XD,N);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FILTER  %%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	m = m0;
	P = P0;
	lh = 0;
    glh = 0;
    dm = 0;
    dP = 0;
	
	
	for k=1:N
        % PREDICTION, p(x_k|y_1:k-1)=N(mean_pred,P_)
        
        % affine transform the unit sigma points

        mean_prev_replicated = repmat(m,1,NS); 
        sig_pred = mean_prev_replicated+chol(P,'lower')*ei;

        
        % propagate through the dynamics function
        sig_pred_propagated=f(sig_pred,k,p);
        % apply the integration formula, wi(1,:) are the weights for mean
        mean_pred = sig_pred_propagated*(wi(1,:)');
        mean_pred_replicated = repmat(mean_pred,1,NS);
        
        if isa(Q,'function_handle')
            QQ = Q(mean_pred,k,p);
        else
            QQ = Q;
        end
        
        % apply the integration formula, wi(2,:) are the weights for covariance
        P_ = (sig_pred_propagated-mean_pred_replicated)*diag(wi(2,:))*...
             (sig_pred_propagated-mean_pred_replicated)'+QQ;
        % compute the cross covariance E[(x_k-1,k-1 - 1-m_k-1,k-1)(x_k,k-1 - m_k,k-1)]
        D = (sig_pred-mean_prev_replicated)*diag(wi(2,:))*...
            (sig_pred_propagated-mean_pred_replicated)';
        
        ms_(:,k) = mean_pred;
        Ps_(:,:,k) = P_;
        Ds(:,:,k) = D;
        
        if k < 5
            PP_ = P_(:);
            PP = P(:);
            fprintf(1,'norm ( %.2f, %.2f, %.2f, %.2f )\n',m'*m,PP'*PP,mean_pred'*mean_pred,PP_'*PP_);
        end
        
        
        % UPDATE, p(x_k|y_1:k)=N(m,P)
        if sum(isnan(y(:,k))) == YD
        	m = mean_pred;
        	P = P_;
        	ms(:,k) = m;
        	Ps(:,:,k) = P;
        	continue;
        end
        % affine transform the unit sigma points
        sig_upd = mean_pred_replicated+chol(P_,'lower')*ei;
        % propagate through the measurement function
        sig_upd_propagated = h(sig_upd,k,p);
        % integrate
        mean_measurement = sig_upd_propagated*wi(1,:)';
        mean_measurement_replicated = repmat(mean_measurement,1,NS);
        
        S = (sig_upd_propagated-mean_measurement_replicated)*diag(wi(2,:))*...
            (sig_upd_propagated-mean_measurement_replicated)'+R;
        C = (sig_upd-mean_pred_replicated)*diag(wi(2,:))*...
            (sig_upd_propagated-mean_measurement_replicated)';
         

        K = C/S;
        d = y(:,k)-mean_measurement;
        m = mean_pred+K*d;
        P = P_-K*S*K';
        ms(:,k) = m;
        Ps(:,:,k) = P;
        
        if no > 0
			lh = lh + likelihood(d,S,YD);
        end
        if no > 1
			% for AR(1) model
            A = p{1};
            H = p{3};
            dA_A = 1;
            dQ_A = 0;
            dR_A = 0;
            
            if k == 1
                m_prev = m0;
                P_prev = P0;
            else
                m_prev = ms(:,k-1);
                P_prev = Ps(:,:,k-1);
            end
            
            dm_ = dA_A*m_prev+A*dm;
            dP_ = dA_A*P_prev*A'+...
                        A*dP*A'+...
                        A*P_prev*dA_A'+...
                        dQ_A;
            % update
            dS = H*dP_*H'+dR_A;
            dK = dP_*H'/S-P_*H'*(S\dS/S);

            dm = dm_+dK*d-K*H*dm_; 
			dP = dP_-dK*S*K'-K*dS*K'-K*S*dK'; 
            
            do0 = trace(S\dS);
            do1 = (H*dm_)'/S*d;
            do2 = d'*(S\dS/S)*d;

            % gradient
            glh = glh + (0.5*(do2-do0)+do1);
            
        end
	end

	if no > 0
		varargout{1} = lh;
	end
	if no > 1
		varargout{2} = glh;
	end
	
end


function lh=likelihood(d,S,YD)
% compute the likelihood addition
    lh = - 0.5*(d'*(S\d) + log(det(S)) + YD*log(2*pi));

end




