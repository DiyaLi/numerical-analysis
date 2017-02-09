function pc =cauchy(g,B,R)
denom = g'*B*g;
T = (norm(g))^3/(R*denom);
if denom <= 0
    T = 1;
else
    T = min(T,1);
end
pc = -T*R*g/norm(g);


