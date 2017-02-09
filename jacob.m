function Jcob = jacob(x01,x02,x03,x04)
syms x1 x2 x3 x4 x5
Jcob = zeros(4,12);
f = x1*sin(x2*x5+ x3) + x4;
Jco = jacobian(f,[x1,x2,x3,x4]);
for m = 1:12
    Jcob(:,m) = subs(Jco,[x1 x2 x3 x4 x5],[x01 x02 x03 x04 m]);
end
Jcob = double(Jcob);