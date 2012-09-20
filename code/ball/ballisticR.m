function [ R ] = ballisticR( lr )

    R = diag(exp(2*lr))*eye(2);

end

