function DotQ = WnbbToDotQ( wnbb,Q )
%使用wnbbs计算四元数的微分
%20170316Salamander
Omega=[ 0,      -wnbb(1),-wnbb(2),  -wnbb(3);
        wnbb(1),0,       wnbb(3),   -wnbb(2);
        wnbb(2),-wnbb(3),0,         wnbb(1);
        wnbb(3),wnbb(2),-wnbb(1),   0];
DotQ=0.5*Omega*Q;



end

