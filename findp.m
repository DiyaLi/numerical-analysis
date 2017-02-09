function p = findp(Q,A,q,b,w0)
p1 = inv([Q A';A zeros(l)])*[-q;b];
l = length(w0);
for k = 1:l
    lamb(k) = p1(2+k);
end
if lamb < 0
    [lamb,mu0] = min(lamb);
    w0(mu0)=[]; A(mu0)=[]; b(mu0)=[];
    p2 = inv([Q A';A zeros(l)])*[-q;b];
end
