function [mm,PP,DD] = SigmaSmoother(mm,PP,mm_,PP_,CC)
	
	XD = size(mm,1);
	N = size(mm,2)-1; % number of measurements
	
	%%%%%%%%%%%%%%%
	% SMOOTHER %%%%
	%%%%%%%%%%%%%%%
	
	
    % The cross-timestep covariance matrices
    DD = zeros(XD,XD,N);
	
    m_s = mm(:,end);
	P_s = PP(:,:,end);
    I2 = 0;    
    
    for k=N:-1:1 % not a typo, last index of mm,PP is N+1
	    % filtering distribution
        m = mm(:,k);
	    P = PP(:,:,k);

	    % prediction distribution
	    m_ = mm_(:,k);
	    P_ = PP_(:,:,k);
	    
	    % the cross-covariance
	    C = CC(:,:,k);	    

	    G = C / P_;
	    
        D = G*P_s;

        m_s = m + G*(m_s-m_);
	    P_s = P +G*(P_s-P_)*G';
	    
                        
	    mm(:,k) = m_s;
	    PP(:,:,k) = P_s;
        DD(:,:,k) = D;
        
        Pj = [PP(:,:,k)       G*PP(:,:,k+1); 
              PP(:,:,k+1)*G'  PP(:,:,k+1)];
        I2 = I2 + Pj + [mm(:,k);mm(:,k+1)]*[mm(:,k);mm(:,k+1)]';
%         try 
%             S = chol(Pj);
%         catch err
%             disp(k);
%             disp(Pj);
%             rethrow(err);
%         end
%                 
%         if k == 2 || k == 3
%             disp('smoother');
%             disp(['between ' num2str(k) ',' num2str(k+1)]);
%             PP(:,:,[k k+1])
%             D
%         end
        
        
    end

    disp(I2);
    
end