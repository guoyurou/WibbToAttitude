function [ Cen,flag ] = LatiLongtiToCen( lati,longti )
%使用载体所在的经纬度求解由地球系e到导航系n的坐标变换矩阵
%20170309 Salamander
%输入单位为弧度
%flag=0表示正常
flag=1;
if lati<-0.5*pi||lati>0.5*pi||longti<-pi||longti>pi
    disp('函数：LatiLongtiToCen 输入经纬度超限')
    return
end

Cen=[-sin(longti),          cos(longti),            0;
    -sin(lati)*cos(longti), -sin(lati)*sin(longti), cos(lati);
    cos(lati)*cos(longti),  cos(lati)*sin(longti),  sin(lati)    ];
flag=0;
end

