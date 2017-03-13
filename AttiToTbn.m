function [Tbn,flag] = AttiToTbn( attitude )
%利用方位俯仰横滚计算由载体系n到导航系b的转换矩阵
%20170309 Salamander
%   attitude=[  psi;
%               theta
%               gama];
% 单位：弧度
%flag=0表示正常
%角度范围检核
Tbn=zeros(3,3);
flag=1;
psi=attitude(1);
theta=attitude(2);
gama=attitude(3);
if psi<0||psi>2*pi||theta<-0.5*pi||theta>0.5*pi||gama<-pi||gama>pi
    disp('函数：AttiToTbn 输入姿态角超限');
    return
else
    Tbn=[cos(gama)*cos(psi)-sin(gama)*sin(theta)*sin(psi),-cos(theta)*sin(psi),sin(gama)*cos(psi)+cos(gama)*sin(theta)*sin(psi);
         cos(gama)*sin(psi)+sin(gama)*sin(theta)*cos(psi), cos(theta)*cos(psi),sin(gama)*sin(psi)-cos(gama)*sin(theta)*cos(psi);
         -sin(gama)*cos(theta),                            sin(theta),         cos(gama)*cos(theta)                            ];
    flag=0;
end
end

