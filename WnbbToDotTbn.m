function DotTbn = WnbbToDotTbn( wnbb,Tbn )
%使用wnbb计算捷联姿态矩阵的微分
%20170319 Salamander

Omega=[ 0,-wnbb(3),wnbb(2);
        wnbb(3),0,-wnbb(1);
        -wnbb(2),wnbb(1),0];
DotTbn=Tbn*Omega;

end

