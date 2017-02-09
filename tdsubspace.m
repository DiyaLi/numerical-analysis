function TS = tdsubspace(g,B,R)
eigval = eig(B);
I= eye(length(eigval));
n = min(eigval);
if n < 0
    a = (-n).*rand(1) - n;
    TS = -inv(B + a*I)*g;
else
    TS = feval(@cauchy,g,B,R);
end