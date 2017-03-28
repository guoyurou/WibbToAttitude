function [Q,flag ] = AttiToQuater( attitude )
%使用姿态角计算四元数
%20170316Salamander
%   attitude=[  psi;
%               theta
%               gama];
% 单位：弧度
%flag=0表示正常
%角度范围检核
Q=zeros(4,1);
flag=1;
psi=attitude(1);
theta=attitude(2);
gama=attitude(3);
if psi<0||psi>2*pi||theta<-0.5*pi||theta>0.5*pi||gama<-pi||gama>pi
    disp('函数：AttiToQuater 输入姿态角超限');
    return
else
    Q=[ cos(psi/2)*cos(theta/2)*cos(gama/2)-sin(psi/2)*sin(theta/2)*sin(gama/2);
        cos(psi/2)*sin(theta/2)*cos(gama/2)-sin(psi/2)*cos(theta/2)*sin(gama/2);
        cos(psi/2)*cos(theta/2)*sin(gama/2)+sin(psi/2)*sin(theta/2)*cos(gama/2);
        cos(psi/2)*sin(theta/2)*sin(gama/2)+sin(psi/2)*cos(theta/2)*cos(gama/2)
        ];
    flag=0;
end
end

