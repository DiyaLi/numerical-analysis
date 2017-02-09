function astar = liness(pk,x0,f0,df0)
a0 = 0; c1 = 10^-4; c2 = 0.9; amax = 10; a = (amax-0).*rand(1) + 0;
theta0 = df0'*pk;
xk = x0;
[fk] = feval(@rf,xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12));
aa = zeros(0); fun = zeros(0); 
aa=cat(1,aa,a0); fun=cat(1,fun,fk);
for i = 2:1000
xk = x0 + a*pk';
[fk,dfk] = feval(@rfdf,xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12));
aa=cat(1,aa,a); fun=cat(1,fun,fk);
    if fk > f0 + c1*a*theta0 || fk >= fun(i-1)
        astar = feval(@zoom,aa(i-1),aa(i),pk,x0,f0,theta0,c1,c2);
        return
    end
    thetak = dfk'*pk;
    if abs(thetak) <= -c2*theta0
        astar = a;
        return
    end
    if thetak >= 0
        astar = feval(@zoom,aa(i),aa(i-1),pk,x0,f0,theta0,c1,c2);
        return
    end
a = (amax-a).*rand(1) + a;
end