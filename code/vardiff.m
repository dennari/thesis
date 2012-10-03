function [ d ] = vardiff( P,S )

  d = 0;
  for k=size(P,3)
    P_ = P(:,:,k);
    S_ = S(:,:,k)*S(:,:,k)';
    d = d + sum((P_(:)-S_(:)).^2);
  end
  d = sqrt(d/size(P,3));


end

