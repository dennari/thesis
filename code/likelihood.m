function lh=likelihood(d,S)
% compute the likelihood addition
    
    lh = d'/S*d + log(det(S)) + size(d,1)*log(2*pi);
    lh = -lh/2;

end

