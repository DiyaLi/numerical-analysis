clear
e = 10^-5; %error 
thetamax = 0.8; 
x = zeros(0); x0 = [-0.45,-0.45,-0.45]; x = cat(1,x,x0); %initial point
fer = zeros(0); fe = 0; 
f0 = feval(@rf,x0); fval = zeros(0); fval = cat(1,fval,f0);
fe = fe + 1;
r = zeros(0); r0 = 1; r = cat(1,r,r0); %initial trust radius
D = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];    
for k=1:1000
   if fval(k,:) < e
       break;
   end
   for i = 1:6
       xnew = x(k,:) + r(k,:)*D(i,:);
       fnew = feval(@rf,xnew); fe = fe + 1;
       rho = 2*(r(k,:))^(3/2);
       if fnew < f0 - rho
            x(k+1,:) = x(k,:) + r(k,:)*D(i,:);
            break;
       end
   end
   x(k+1,:) = x(k,:) + r(k,:)*D(i,:);   
   fnew = feval(@rf,x(k+1,:)); fe = fe + 1;
   rho = 2*(r(k,:))^(3/2);
   if fnew < f0 - rho     
       fval = cat(1,fval,fnew);     
       ph = 1 + rand(1);
       r = cat(1,r,r(k,:)*ph);
   else
       x(k+1,:) = x(k,:);
       fval = cat(1,fval,fval(k)); 
       theta = thetamax*rand(1);
       r = cat(1,r,r(k,:)*theta);
   end
   fer = cat(1,fer,fe);
   fe = 0;
end
iteration = linspace(1,length(r),length(r));
fer = cat(1,fer,0);
figure
subplot(2,1,1)   
plot(iteration,fval-1)
title('The error at each step versus the number of iterations')
xlabel('Iterations')
ylabel('Error')
subplot(2,1,2) 
plot(iteration,fer)
title('The number of function evaluations versus the number of iterations')
xlabel('Iterations')
ylabel('Number of function evaluations')


