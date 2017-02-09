function a = backss(pk,x0,f0,df0)
c = 10^-4; rho = 0.9 ; a = 1;
xk = x0 + a*pk';
fk = feval(@rf,xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12));
while fk > f0 + c*a*df0'*pk
    a = a*rho;
    if a < 0.0001
        return
    end
    xk = x0 + a*pk';
    fk = feval(@rf,xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12));
end
