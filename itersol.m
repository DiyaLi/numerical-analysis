function IS = itersol(g,B,R)
eigval = eig(B);
I = eye(length(eigval));
n = -min(eigval)+0.3;
p = -inv(B + n*I)*g; 
e = 1e-8;
while norm(p)-R > e
       C = chol(B + n*I);
       pl = (C'*C)\(-g);
       ql = C'\pl;
       n = n + (norm(pl)/norm(ql))^2*(norm(pl)-R)/R;
       p = -inv(B + n*I)*g;
end
IS = -inv(B+n*I)*g;