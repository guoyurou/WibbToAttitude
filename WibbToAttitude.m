%载体不做位移，通过陀螺仪测量值解算姿态角
%20170309 Salamander
%程序计算中所有角度一律使用弧度为单位，只在输入输出阶段涉及度与弧度转换
%20170316 引入四元数处理,使用fscanf按行读入数据而不是一次性全部读入，模拟实机运行情况
%20170322 开始处理经纬度与东向北向速度
%20170329 开始处理高度通道
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
vn=[0;0;0];
vu=0;%高度通道单独初始化
g=AccuG(lati,h);


fileID = fopen('E:\Documents\惯导实验\data20110913\imu_ENU.txt','r');
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';
[IMUSource,FlagRead]=fscanf(fileID,formatSpec);
if FlagRead~=7
    disp('数据读取失败')
    stop
end
wibb=[IMUSource(5);IMUSource(6);IMUSource(7)];
nepoch=65000;
ResOffer=zeros(nepoch,10);
ResOffer=mat2dataset(ResOffer,'VarNames',{'time','psi','theta','gama','ve','vn','vu','lati','longti','h'});

DeltaT=0.01;
DeltaTOfTbn=2;%Tbn矩阵即姿态的更新周期=DeltaTOfTbn*DeltaT 使用四阶龙格库塔法，所以此处姿态数据的更新周期是陀螺仪输出周期的2倍
DeltaTOfV=10;%每间隔DeltaTOfV个采样间隔更新一次东向北向速度
DeltaTOfNoQ=10;%每隔DeltaTOfNoQ个采样间隔归一化一次四元数
DeltaTOfPo=50;%每隔DeltaTOfPo个采样间隔解算一次位置
DeltaTOfH=50;%每隔DeltaTOfH个采样间隔更新一次高度与天向速度
FlagMiddleRefresh=0;
K_two=0.01;%高度通道二阶阻尼回路系数
K_one=1;%高度通道二阶阻尼回路系数
high_c=importHighfile('E:\Documents\惯导实验\data20110913\high.txt');%外部高度信息用于实现高度阻尼

%% 计算初始winb wnbb Tbn

%求wenb
Rm=Re*(1-2*f+3*f*(sin(lati))^2);
Rn=Re*(1+f*(sin(lati))^2);
wenn=[  -vn(2)/(Rm+h);
        vn(1)/(Rm+h);
        vn(1)*tan(lati)/(Rn+h)];
[Tbn,FlagFunc]=AttiToTbn(AttitudeIni);%由b到n的转换矩阵，b为下标n为上标
[Q,FlagFun]=AttiToQuater(AttitudeIni);%由姿态角计算四元数Q
if FlagFunc~=0
    stop
end
wenb=Tbn'*wenn;
%求wieb
Cen=LatiLongtiToCen(lati,longti);
DotCen=WennToDotCen(wenn,Cen);
wien=Cen*wiee;
wieb=Tbn'*wien;
%求winb
winb=wieb+wenb;
wnbb_pre=wibb-winb;
%求初始导航系下的加速度DotVn_pre
fb_ini=[IMUSource(2);IMUSource(3);IMUSource(4)];
fn_ini=Tbn*fb_ini;
DotVn=fn_ini-cross(2*wien+wenn,vn)+[0;0;-g];%此处g之前是否应当加负号？
DotVu=DotVn(3);%高度通道单独初始化
%% 开始读取观测
for i=1:nepoch
    [IMUSource,FlagRead]=fscanf(fileID,formatSpec);
    if FlagRead~=7
        info=sprintf('%s%d','数据读取失败或读到文件尾，结束循环，已读取行数：',i);
        disp(info)
        fclose(fileID);
        break
    end
    wibb=[IMUSource(5);IMUSource(6);IMUSource(7)];
    fb=[IMUSource(2);IMUSource(3);IMUSource(4)];
    ResOffer.time(i)=IMUSource(1);
    %% 更新姿态
    
    if mod(ResOffer.time(i),DeltaTOfTbn)~=0
        wnbb_middle=wibb-winb;
        FlagMiddleRefresh=1;
    elseif FlagMiddleRefresh==1
        wnbb_now=wibb-winb;
        %四阶龙哥库塔求四元数Q
        K1=WnbbToDotQ(wnbb_pre,Q);
        K2=WnbbToDotQ(wnbb_middle,Q+K1*DeltaT*DeltaTOfTbn/2);
        K3=WnbbToDotQ(wnbb_middle,Q+K2*DeltaT*DeltaTOfTbn/2);
        K4=WnbbToDotQ(wnbb_now,Q+K3*DeltaT*DeltaTOfTbn);
        Q=Q+DeltaT*DeltaTOfTbn*(K1+2*K2+2*K3+K4)/6;
        %四元数归一化
        if mod(ResOffer.time(i),DeltaTOfNoQ)==0
            Q=Q/norm(Q);
        end

        %由四元数计算Tbn
        Tbn=QToTbn(Q);
        %输出姿态角
        [ResOffer.psi(i),ResOffer.theta(i),ResOffer.gama(i),FlagFunc]=TbnToAttitude(Tbn);
        if FlagFunc~=0
            disp('Tbn矩阵不可靠，姿态角未更新')  
            if(i<=DeltaTOfTbn)
                ResOffer(i,2:4)=AttitudeIni;
            else
                ResOffer(i,2:4)=ResOffer(i-DeltaTOfTbn,2:4);
            end
        end
        fn=Tbn*fb;
    end
    %% 更新速度
    
    if mod(ResOffer.time(i),DeltaTOfV)==0


        %二阶龙格库塔法更新速度
        K1=DotVn;
        wenn=[  -(vn(2)+K1(2)*DeltaTOfV*DeltaT)/(Rm+h);
        (vn(1)+K1(1)*DeltaTOfV*DeltaT)/(Rm+h);
        (vn(1)+K1(1)*DeltaTOfV*DeltaT)*tan(lati)/(Rn+h)];
        K2=fn-cross(2*wien+wenn,vn+K1*DeltaTOfV*DeltaT)+[0;0;-g];
        vn=vn+(K1+K2)*DeltaTOfV*DeltaT/2;
        %更新数据
        wenn=[  -vn(2)/(Rm+h);
                vn(1)/(Rm+h);
                vn(1)*tan(lati)/(Rn+h)];
        DotVn=fn-cross(2*wien+wenn,vn)+[0;0;-g];
        ResOffer.ve(i)=vn(1);
        ResOffer.vn(i)=vn(2);

    end
    %% 更新经纬度
    
    if mod(ResOffer.time(i),DeltaTOfPo)==0
        Cen=Cen+DeltaTOfPo*DeltaT*DotCen;
        DotCen=WennToDotCen(wenn,Cen);
        [lati,longti]=CenToLaLon(Cen);
        %更新其他数据
        Rm=Re*(1-2*f+3*f*(sin(lati))^2);
        Rn=Re*(1+f*(sin(lati))^2);
        wenn=[  -vn(2)/(Rm+h);
        vn(1)/(Rm+h);
        vn(1)*tan(lati)/(Rn+h)];
        wien=Cen*wiee;
        ResOffer.lati(i)=lati;
        ResOffer.longti(i)=longti;
        
    end
        
    %% 更新wnbb_pre
    if mod(ResOffer.time(i),DeltaTOfTbn)==0&&FlagMiddleRefresh==1
        winb=Tbn'*(wien+wenn);
        wnbb_pre=wibb-winb;
        FlagMiddleRefresh=0;
    end
    
    %% 更新高度通道
    if mod(ResOffer.time(i),DeltaTOfH)==0
        index_hc=find(high_c.time==ResOffer.time(i));
        if size(index_hc)==1
            hc=high_c.h(index_hc);
            %一阶欧拉法更新高度
            h=h+(vu-K_one*(h-hc))*DeltaTOfH*DeltaT;
            ResOffer.h(i)=h;
            vu=vu+(DotVu-K_two*(h-hc))*DeltaTOfH*DeltaT;
            ResOffer.vu(i)=vu;
            DotVu=fn(3)+(2*wien(2)+wenn(2))*vn(1)-(2*wien(1)+wenn(1))*vn(2)-AccuG(lati,h);
            vn(3)=vu;
            DotVn(3)=DotVu;
            %此处若size(index_hc)==1不成立，则下次数值积分时其实步长就不再是DeltaTOfH*DeltaT，应当累加，此处暂且不考虑该情况
        end
    end
end
%输出角度以度为单位    
ResOffer.psi=ResOffer.psi*180/pi;
ResOffer.theta=ResOffer.theta*180/pi;
ResOffer.gama=ResOffer.gama*180/pi;
ResOffer.lati=ResOffer.lati*180/pi;
ResOffer.longti=ResOffer.longti*180/pi;
save('ResOfWibbToAttitude','ResOffer');

%% 比对程序WibbToAttitude计算得到的结果与ins_c.txt中的参考姿态数据
INSSource=imporINSCtfile('E:\Documents\惯导实验\data20110913\ins_c.txt');
index=(ResOffer.psi~=0);
subplot(3,3,1)
plot(ResOffer.time(index),ResOffer.psi(index),'.');
hold on
plot(INSSource.time,INSSource.psi,'.r')
hold off
subplot(3,3,2)
plot(ResOffer.time(index),ResOffer.theta(index),'.');    
hold on
plot(INSSource.time,INSSource.theta,'.r')
hold off
subplot(3,3,3)
plot(ResOffer.time(index),ResOffer.gama(index),'.');
hold on
plot(INSSource.time,INSSource.gama,'.r')
hold off

index_v=(ResOffer.ve~=0);
subplot(3,3,4)
plot(ResOffer.time(index_v),ResOffer.ve(index_v),'.');
hold on
plot(INSSource.time,INSSource.Ve,'.r')
hold off
subplot(3,3,5)
plot(ResOffer.time(index_v),ResOffer.vn(index_v),'.');
hold on
plot(INSSource.time,INSSource.Vn,'.r')
hold off
index_vu=(ResOffer.vu~=0);
subplot(3,3,6)
plot(ResOffer.time(index_vu),ResOffer.vu(index_vu),'.');
hold on
plot(INSSource.time,INSSource.Vu,'.r')
hold off

index_p=(ResOffer.lati~=0);
subplot(3,3,7)
plot(ResOffer.time(index_p),ResOffer.lati(index_p),'.');
hold on
plot(INSSource.time,INSSource.lati,'.r')
hold off
subplot(3,3,8)
plot(ResOffer.time(index_p),ResOffer.longti(index_p),'.');
hold on
plot(INSSource.time,INSSource.longti,'.r')
hold off
index_h=(ResOffer.h~=0);
subplot(3,3,9)
plot(ResOffer.time(index_h),ResOffer.h(index_h),'.');
hold on
plot(INSSource.time,INSSource.high,'.r')
hold off


