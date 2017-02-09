clear
i = 1; Rh = 1; 
R =0.849129305868777; mu = 0.008927919643547;
fmin = 3.62059518154553*(10^-19); e = 10^-8;
error = zeros(0); iter = zeros(0); f = zeros(0); cput = zeros(0);
x = [1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2 1.2];
[fx,df,ddf] = feval(@rnumf,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12));
f = cat(1,f,fx);
while norm(df) > e
   t = cputime;
%   pk = feval(@dogleg,df,ddf,R);
%   pk = feval(@tdsubspace,df,ddf,R);
   pk = feval(@itersol,df,ddf,R);
   xk = x + pk';
   fk = feval(@rf,xk(1),xk(2),xk(3),xk(4),xk(5),xk(6),xk(7),xk(8),xk(9),xk(10),xk(11),xk(12)); 
   m = -df'*pk - pk'*ddf*pk/2;
   rho = (f(i) - fk)/m;
   if rho < 1/4
       R = R/4;
   else
       if (rho > 3/4) && (norm(pk) == R)
           R = min(2*R,Rh);
       end
   end
   if rho > mu
       x = x + pk';
       [fx,df,ddf] = feval(@rnumf,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12));
       f = cat(1,f,fx);
   else
       f = cat(1,f,fx);
   end
   iter = cat(1,iter,i);
   error = cat(1,error,abs(f(i+1) - fmin));   
   cput = cat(1,cput,cputime-t); 
   i = i + 1;
end
figure
subplot(2,1,1)   
semilogy(iter,error)
title('The logarithm of the error at each step versus the number of iterations')
xlabel('Iterations')
ylabel('Log(Error)')
subplot(2,1,2) 
semilogy(cput,error,'.')
title('The logarithm of the error at each step versus the CPU time')
xlabel('CPU time (s)')
ylabel('Log(Error)')

