function [ Q,A ] = myQR( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = size(A,2);
m = size(A,1);
V = zeros(m,n);
Q = zeros(m,m);
for k=1:n
    x = A(k:m,k);
    s = sign(x(1)) || 1;
    e = zeros(size(x)); e(1) = 1;
    v = s*norm(x)*e+x;
    v = v/norm(v);
    A(k:m,k:n) = (eye(m-k+1)-2*(v*v'))*A(k:m,k:n);
    V(k:m,k) = v;
end
e = zeros(m,1);
for i=1:m
    e(i) = 1;
    for k=n:-1:1
        v = V(k:m,k);
        e(k:m) = e(k:m)-2*v*(v'*e(k:m));
    end
    Q(:,i) = e;
    e = 0*e;
end
