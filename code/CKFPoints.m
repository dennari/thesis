function [ei,wi] = CKFPoints(XD)

	% CKF unit sigma points and weights
	ei = sqrt(XD)*[eye(XD) -eye(XD)];
	wi = ones(1,2*XD);
	wi = wi/sum(wi);
	
		
end