function DL = dogleg(g,B,R)
pB = -inv(B)*g;
if R >= norm(pB)
    DL = pB;
else
    pu = -g'*g/(g'*B*g)*g;
    syms t
    equ = norm(pu+(t-1)*(pB-pu))==R;
    t = double(solve(equ,t));
    T = double(t(t>0));
    T = max(T);
    if 0 <= T && T <= 1
        DL = T*pu;
    else if 1 <= T && T <= 2
            DL = pu + (T - 1)*(pB - pu);
         end
    end    
end
