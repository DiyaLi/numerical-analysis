clear
%Model 1
%val = [42 45 53 63 71 79 82 81 75 65 56 44]; x0 = [5 0 5 50];
%Model 2
%val = [21 24 30 38 47 56 60 59 51 39 31 24];x0 = [5 0 5 30];
%Model 3
val = [3.07 2.8 3.62 3.46 4.33 4.02 4.25 3.58 3.11 2.8 2.87 2.95];
x0 = [1 1 1 1];
x = zeros(0);  x = cat(1,x,x0); f = zeros(0); r = zeros(0);
iter = zeros(0); iter0 = 0; err = zeros(0); cput = zeros(0);
for m = 1:12
    f0(m) = feval(@model, x(1,1), x(1,2), x(1,3), x(1,4),m);
    r0(m) = f0(m) - val(m);
end
J = feval(@jacob, x0(1), x0(2), x0(3), x0(4)); lamb = 0.001;
while norm(J*r0') > 10^(-8) 
    cpu = cputime;
    g = J*r0';
    H = J*J'+lamb*eye(4);
    p = -inv(H)*g;
    xk = x0 + p';
    for m = 1:12
        fk(m) = feval(@model, xk(1), xk(2), xk(3), xk(4),m);
        rk(m) = fk(m) - val(m);
    end
    if norm(rk) < norm(r0)
        x0 = xk; f0 = fk; r0 = rk; lamb = lamb/10;
    else
        x0 = x0; f0 = f0; r0 = r0; lamb = lamb*10;
    end
    iter0 = iter0 + 1;
    err0 = norm(r0);
    x = cat(1,x,x0); f = cat(1,f,f0); r = cat(1,r,r0); 
    iter = cat(1,iter,iter0); err = cat(1,err,err0); 
    J = feval(@jacob, x0(1), x0(2), x0(3), x0(4));
    cput = cat(1,cput,cputime-cpu); 
end 
m = 1:12;
t = linspace(1,12,1000);
y = x0(1)*sin(x0(2)*t + x0(3)) + x0(4);
subplot(3,1,1) 
plot(t,y,'-',m,val,'o')
xlabel('Month')
ylabel('Av. Precipitation.')
subplot(3,1,2) 
semilogy(iter,err)
xlabel('Iteration')
ylabel('Log(Error)')
subplot(3,1,3) 
semilogy(cput,err)
xlabel('Time of CPU')
ylabel('Log(Error)')

