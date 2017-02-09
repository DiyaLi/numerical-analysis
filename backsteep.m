clear
a0 = 1; %step length
c = 10^-4; rho = 0.8; %backtrack step
a = a0; i = 1; e = 10^-8;
itermax = 10000;
x0 = [1.2 1.2]; % initial point
%x0 = [-1.2 1]; 
xsize = length(x0);
[f0,df0,ddf0] = feval(@rnumf,x0(1),x0(2));
storex = zeros(itermax,xsize); % record x value
storef = zeros(itermax,1); % record f value
storendf = zeros(itermax,1); % record norm df value
storess = zeros(itermax,1); % record step size 
storex(i,:)=x0; storef(i,:)=f0; storendf(i,:)=norm(df0); storess(i,:)=a;
while norm(df0) > e
    if i < 100
        pk = feval(@steepdir,xsize,df0);
        % pk = feval(@numdir,df0,ddf0); 
        a = a0;
        a = feval(@backss,pk,x0,c,a,f0,df0,rho);
        xk = x0 + a*pk'; [fk,dfk,ddfk] = feval(@rnumf,xk(1),xk(2));
        x0 = xk; f0 = fk; df0 = dfk; ddf0 = ddfk;
        i = i + 1;
        ndf = norm(df0);
        storex(i,:)=x0; 
        storef(i,:)=f0; 
        storendf(i,:)=ndf; 
        storess(i,:)=a;
    else 
        return
    end
end
minimizef = storef(i,:);
error = zeros(i,1); % record errors
func = zeros(i,1); % record funtion value
iter = zeros(i,1);
for j=1:i
    iter(j,1) = j;
    func(j,1) = storef(j,:);
    error(j,1)= abs(minimizef - func(j,1));
end
plot(iter,error)

    
