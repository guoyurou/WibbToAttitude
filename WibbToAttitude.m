%载体不做位移，通过陀螺仪测量值解算姿态角
%20170309 Salamander
%程序计算中所有角度一律使用弧度为单位，只在输入输出阶段涉及度与弧度转换

clc
clear

%% 预置固定参数
Re=6378137.0;
f=1.0/298.257;
wiee=[  0;
        0;
        7.2722e-5];

%% 初始值
AttitudeIni=[89.850004;1.9113951;1.0572407]*pi/180;
lati=39.813330*pi/180;
longti=116.15326*pi/180;
h=70.0;
vn=[0,0,0];

IMUSource=importIMUfile('E:\Documents\惯导实验\导航系统设计程序数据\data20110913\imu_ENU.txt');
wibb_ini=[-1.7630406e-003;-1.2578430e-004;-2.6035600e-005];
nepoch=size(IMUSource);
AttitudeOffer=zeros(ceil(nepoch(1)/2),4);
DeltaT=0.02;%由于使用四阶龙格库塔法，所以此处的数据输出周期DeltaT其实是陀螺仪输出周期的2倍
FlagMiddleRefresh=0;

%%参数
method=4;%微分方程解法选项
%% 计算初始winb wnbb Tbn

%求wenb
Rm=Re*(1-2*f+3*f*(sin(lati))^2);
Rn=Re*(1+f*(sin(lati))^2);
wenn=[  -vn(2)/(Rm+h);
        vn(1)/(Rm+h);
        vn(1)*tan(lati)/(Rn+h)];
[Tbn_pre,FlagFunc]=AttiToTbn(AttitudeIni);%由b到n的转换矩阵，b为下标n为上标
if FlagFunc~=0
    stop
end
wenb=Tbn_pre'*wenn;
%求wieb
Cen=LatiLongtiToCen(lati,longti);
wien=Cen*wiee;
wieb=Tbn_pre'*wien;
%求winb
winb_pre=wieb+wenb;
wnbb_pre=wibb_ini-winb_pre;

%% 开始读取观测
for i=1:nepoch(1)
    wibb=[IMUSource.wibbx(i);IMUSource.wibby(i);IMUSource.wibbz(i)];
    if mod(i,2)~=0
        wnbb_middle=wibb-winb_pre;
        FlagMiddleRefresh=1;
        continue
    end
    AttitudeOffer(i/2,1)=IMUSource.time(i);
    wnbb_now=wibb-winb_pre;
    %四阶龙哥库塔求Tbn
    if FlagMiddleRefresh==0
        disp('四阶龙格库塔法无法实施，wnbb_middle值未更新')
        continue
    end
    if method==4
        
        K1=WnbbToDotTbn(wnbb_pre,Tbn_pre);
        K2=WnbbToDotTbn(wnbb_middle,Tbn_pre+K1*DeltaT/2);
        K3=WnbbToDotTbn(wnbb_middle,Tbn_pre+K2*DeltaT/2);
        K4=WnbbToDotTbn(wnbb_now,Tbn_pre+K3*DeltaT);
        Tbn_now=Tbn_pre+DeltaT*(K1+2*K2+2*K3+K4)/6;
    end
   
    %更新下一次循环使用的数据,此处静态数据暂不考虑wien的更新
    winb_pre=Tbn_now'*wien+Tbn_now'*wenn;
    wnbb_pre=wibb-winb_pre;
    Tbn_pre=Tbn_now;
    FlagMiddleRefresh=0;

    %输出姿态角
    [AttitudeOffer(i/2,2:4), FlagFunc]=TbnToAttitude(Tbn_now);
    if FlagFunc~=0
        disp('Tbn矩阵不可靠，姿态角未更新')  
        if(i<=2)
            AttitudeOffer(i/2)=AttitudeIni;
        else
            AttitudeOffer(i/2)=AttitudeOffer(i/2-1);
        end
     end
    
    
end
%输出角度以度为单位    
AttitudeOffer(:,2:4)=AttitudeOffer(:,2:4)*180/pi;
save('ResOfWibbToAttitude','AttitudeOffer');

%% 比对程序WibbToAttitude计算得到的姿态与ins_c.txt中的参考姿态数据


INSSource=imporINSCtfile('E:\Documents\惯导实验\导航系统设计程序数据\data20110913\ins_c.txt');
index=(AttitudeOffer(:,1)~=0);
subplot(3,1,1)
plot(AttitudeOffer(index,1),AttitudeOffer(index,2),'.');
hold on
plot(INSSource.time,INSSource.psi,'.r')
hold off
subplot(3,1,2)
plot(AttitudeOffer(index,1),AttitudeOffer(index,3),'.');    
hold on
plot(INSSource.time,INSSource.theta,'.r')
hold off
subplot(3,1,3)
plot(AttitudeOffer(index,1),AttitudeOffer(index,4),'.');
hold on
plot(INSSource.time,INSSource.gama,'.r')
hold off

