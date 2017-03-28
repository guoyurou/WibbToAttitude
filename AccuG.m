function g = AccuG( lati,h )
%使用1930年国际重力公式计算重力加速度g
%20170322Salamander
%输入：lati 纬度 单位 弧度
%       h   高 单位 米
%输出单位：米/s^2
g=9.78049*(1+0.005288*sin(lati)^2-0.000006*sin(2*lati)^2-0.0003086*h);


end

