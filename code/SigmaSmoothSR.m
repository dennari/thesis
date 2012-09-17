function [ms,Ss,G] = SigmaSmoothSR(ms,Ss,f,SQ,u,usig,w)
	
	XD = size(ms,1);
	N = size(ms,2);
  NS = size(usig,2);
  G = zeros(XD,XD,N);
  if nargin < 5 || isempty(u)
    u = zeros(XD,1);
  end
  if size(w,2) < 2
    w(:,2) = w(:,1);
  end
  wm = w(:,1);
  wp = w(:,2);
  
  for k=N-1:-1:1

    
    sig = Ss(:,:,k)*usig;
    CW1 = sig*diag(wp);

    sig = f(sig+repmat(ms(:,k),1,NS))+repmat(u,1,NS);
    m_ = sig*wm;
    CW2 = (sig-repmat(m_,1,NS))*diag(wp);


    [~,U] = qr([CW2 SQ; CW1 zeros(XD,XD)]',0);
    U = U';
    
    
    % smoother gain
    G(:,:,k) = U(XD+1:end,1:XD)/U(1:XD,1:XD);
    ms(:,k) = ms(:,k) + G(:,:,k)*(ms(:,k+1)-m_);

    [~,S] = qr([U(XD+1:end,XD+1:end) G(:,:,k)*Ss(:,:,k+1)]',0);
    Ss(:,:,k) = S';


  end


    
end