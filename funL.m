function L = funL(x1,x2,x3,x4,x5,b1,b2,b3)
syms a b c d e m1 m2 m3
f = exp(a*b*c*d*e)-0.5*(a^3+b^3+1)^2;
df = jacobian(f,[a,b,c,d,e]);
c = [a^2+b^2+c^2+d^2+e^2-10; b*c-5*d*e;a^3+b^3+1];
dc = jacobian(c,[a,b,c,d,e]);
A = df-dc'*[m1; m2; m3];
dA = jacobian(A,[a,b,c,d,e]);
L = subs(L,[a,b,c,d,e,m1,m2,m3],[x1,x2,x3,x4,x5,b1,b2,b3]);
L = double(L);