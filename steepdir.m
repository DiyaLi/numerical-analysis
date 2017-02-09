function pk = steepdir(xsize,df0)
Bk = eye(xsize);
pk = -inv(Bk)*df0;