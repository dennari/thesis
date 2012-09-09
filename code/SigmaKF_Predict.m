function [m,P,C] = SigmaKF_Predict(m,P,f,Q,usig,w)



	% the number of sigma points
	NS = size(usig,2);

	
    % PREDICTION, p(x_k|y_1:k-1)=N(mean_pred,P_)


    m_rep = repmat(m,1,NS); 

    sig = m_rep+chol(P,'lower')*usig;


    % propagate through the dynamics function
    sigp=f(sig);
    % apply the integration formula, wi(1,:) are the weights for mean
    m = sigp*w;
    m_pred_rep = repmat(m,1,NS);

    % center
    d = sigp-m_pred_rep;
    % weight
    P = d*diag(w)*d'+Q;
    % compute the cross covariance E[(x_k-1,k-1 - 1-m_k-1,k-1)(x_k,k-1 - m_k,k-1)]
    C = (sig-m_rep)*diag(w)*(sigp-m_pred_rep)';
                  
end




