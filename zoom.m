function astar = zoom(alow,ahigh,pk,x0,f0,theta0,c1,c2)
for i=1:10000
    aj = (alow + ahigh)/2;
    xk = x0 + aj*pk';
    xlow = x0 + alow*pk';
    [fk,dfk] = feval(@rfdf,xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12));
    [flow] = feval(@rf,xlow(1),xlow(2),xlow(3),xlow(4),xlow(5),xlow(6),xlow(7),xlow(8),xlow(9),xlow(10),xlow(11),xlow(12));
    if fk > f0 + c1*aj*theta0 || fk >= flow 
        ahigh = aj;
    else
        thetak = dfk'*pk;
        if abs(thetak) <= -c2*theta0
            astar = aj; 
        return
        end
        if thetak*(ahigh - alow) >= 0
            ahigh = alow;
        end
        alow = aj;
    end
end

