function [ lati,lonti ] = CenToLaLon( Cen )
%使用由地球系e转到导航系n的方向余弦矩阵Cen解算纬度lati经度lati，指北方位系统，不考虑自由方位角
%Cen指北方位系统方位余弦矩阵
%lati 纬度 单位：弧度
%lonti 经度 单位：弧度
%后续补充对Cen正交性的检验
lati=asin(Cen(3,3));
lonti=atan(Cen(3,2)/Cen(3,1));
if Cen(3,1)<0
    if lonti>0
        lonti=lonti-pi;
    else
        lonti=lonti+pi;
    end
end



end

