function [m,P,S,d] = SigmaKF_Update(m,P,y,h,R,usig,w)


    
    % the number of sigma points
    NS = size(usig,2);

    % UPDATE, p(x_k|y_1:k)=N(m,P)

    m_rep = repmat(m,1,NS); 
    sig = m_rep+chol(P,'lower')*usig;


    % propagate through the measurement function
    sigp=h(sig);
    % apply the integration formula, wi(1,:) are the weights for mean
    ym = sigp*w;
    ym_rep = repmat(ym,1,NS);

    d = sigp - ym_rep;
    S = d*diag(w)*d'+R;
    C = (sig-m_rep)*diag(w)*(sigp-ym_rep)';

    % Kalman gain
    K = C/S;
    % residual
    d = y-ym;
    % updated mean
    m = m+K*d;
    % updated covariance matrix
    P = P-K*S*K';
      
	





