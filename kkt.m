function [d , lamb, mu] = kkt(Q,A,B,g)
[c,k] = size(A); lamb = zeros(0); mu = zeros(0);
M = [Q A' B';A zeros(c) zeros(c);B zeros(c) zeros(c)];
m = [-g; zeros(c,1); zeros(c,1)];
v = inv(M)*m;
d = [v(1);v(2)];
for u=1:c
    lamb = cat(1,lamb,v(2+u));
    mu = cat(1,mu,v(2+c+u));
end