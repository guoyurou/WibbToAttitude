function DotCen = WennToDotCen( wenn,Cen )
%使用wenn与Cen计算Cen的变化率
%20170323SalamandeR
DotCen=-[0,-wenn(3),wenn(2);
        wenn(3),0,-wenn(1);
        -wenn(2),wenn(1),0]*Cen;
end

