function [ psi,theta,gama,flag ] = TbnToAttitude( Tbn )
%使用由载体系b转化为导航系n所需的Tbn矩阵计算姿态角
%20170310Salamander
%flag=0表示正常

%后续补充对输入的Tbn的正交性实施检验
psi=atan(-Tbn(1,2)/Tbn(2,2));
if Tbn(2,2)>0
    if psi<0
        psi=psi+2*pi;
    end
else
    psi=psi+pi;
end
theta=asin(Tbn(3,2));
gama=atan(-Tbn(3,1)/Tbn(3,3));

flag=0;


end

