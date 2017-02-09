clear
x0 = [-1.8, 1.7, 1.9, -0.8, -0.8];
xstar = [-1.71, 1.59, 1.82,-0.763,-0.763];
lamb = [1 1 1];
dfk = 1;
for k = 1:20
    fk = feval(@fun, x0(1), x0(2), x0(3), x0(4), x0(5));
    dfk = feval(@dfun, x0(1), x0(2), x0(3), x0(4), x0(5));
    dL = feval(@funL, x0(1), x0(2), x0(3), x0(4), x0(5),lamb(1),lamb(2),lamb(3));
    c = feval(@funcc, x0(1), x0(2), x0(3), x0(4), x0(5));
    A = feval(@func, x0(1), x0(2), x0(3), x0(4), x0(5));

    N = [dL -A';A zeros(3)];
    n = [-dfk'; -c];
    r = N\n;
    for m=1:5
        p(m)=r(m);
    end
    for m=6:8
        lm(m-5)=r(m);
    end
    x0 = x0 + p;
    lamb = lm;
    k = k + 1;
end


