clear all
close all
clc
%% 仿真参数设置

%----Global constants----%
global N m l j cn ct  K_P K_D K_theta alpha  omega delta  H
%-----仿真时间、步长 参数设置----%
t0=0; tfin=60; tsamp=0.1; num_seg_t=(tfin-t0)/tsamp; %%起始时间：0，结束时间：45，取样时间：0.1，取样数：450个

%-----蛇形机器人参数设置：长度、质量----%
N=10;m=1; l=0.14; j=0.0016; 

%-----地面环境参数设置：法相、切向摩擦系数----%
ct=1; cn=10; %Links' parameters

%-----关节角控制器参数设置：PD控制器----%
K_D=5; K_P=20;K_theta=0.30;

%-----关节摆动函数参数设置----%
alpha=20*pi/180;omega=pi/2;delta=2*pi/(N-1);

sum_py=0;
%%

 phi1=zeros(num_seg_t,1); phi2=zeros(num_seg_t,1); phi3=zeros(num_seg_t,1); phi4=zeros(num_seg_t,1); phi5=zeros(num_seg_t,1);
 phi6=zeros(num_seg_t,1); phi7=zeros(num_seg_t,1); phi8=zeros(num_seg_t,1); phi9=zeros(num_seg_t,1);
 
 theta10=zeros(num_seg_t,1);
 px=zeros(num_seg_t,1);px=zeros(num_seg_t,1);
 
 d_phi1=zeros(num_seg_t,1);d_phi2=zeros(num_seg_t,1);d_phi3=zeros(num_seg_t,1);d_phi4=zeros(num_seg_t,1);d_phi5=zeros(num_seg_t,1);
 d_phi6=zeros(num_seg_t,1);d_phi7=zeros(num_seg_t,1);d_phi8=zeros(num_seg_t,1);d_phi9=zeros(num_seg_t,1);
 
 d_theta10=zeros(num_seg_t,1);
 d_px=zeros(num_seg_t,1);d_px=zeros(num_seg_t,1);


init_config_vars=[(10/180)*pi*ones(9,1);zeros(13,1);0;0.1]; %initial condition of configuration vars 初始状态向量;起始位置（0，3）

H=[1,1,1,1,1,1,1,1,1,1
   0,1,1,1,1,1,1,1,1,1
   0,0,1,1,1,1,1,1,1,1
   0,0,0,1,1,1,1,1,1,1
   0,0,0,0,1,1,1,1,1,1
   0,0,0,0,0,1,1,1,1,1
   0,0,0,0,0,0,1,1,1,1
   0,0,0,0,0,0,0,1,1,1
   0,0,0,0,0,0,0,0,1,1
   0,0,0,0,0,0,0,0,0,1];

global  last_phi last_dphi u_bar 

%initialize vector
ref_phi=zeros(N-1,1); dref_phi=zeros(N-1,1); ddref_phi=zeros(N-1,1); u_bar=zeros(N-1,1);
P_y=zeros(num_seg_t,1);
for k=1:1:num_seg_t
    k
% time interval 时间间隔
    current_t=t0+k*tsamp;
    interval=[current_t-tsamp,current_t]; %[already known, haven't known] the output of the system  
%last output
last_phi=[init_config_vars(1);init_config_vars(2);init_config_vars(3);init_config_vars(4);init_config_vars(5);init_config_vars(6);init_config_vars(7);init_config_vars(8);init_config_vars(9)];
last_dphi=[init_config_vars(13);init_config_vars(14);init_config_vars(15);init_config_vars(16);init_config_vars(17);init_config_vars(18);init_config_vars(19);init_config_vars(20);init_config_vars(21)];
last_px=init_config_vars(11);
last_py=init_config_vars(12);

last_theta=H*[last_phi;init_config_vars(10)];
sum_last_theta=last_theta(1)+last_theta(2)+last_theta(3)+last_theta(4)+last_theta(5)+last_theta(6)+last_theta(7)+last_theta(8)+last_theta(9)+last_theta(10);
last_theta_average=sum_last_theta/N;
% y_ref=5;

Delta=0.05; %影响跟踪轨迹的速度


theta_ref_bar=-atan((last_py)/Delta);
delta_theta=rem((last_theta_average-theta_ref_bar),pi);
phi0=K_theta*(delta_theta);

%%给定关节运动参考角度
    for i=1:(N-1)  
       ref_phi(i)=alpha*sin(omega*current_t+i*delta)+phi0;
       dref_phi(i)=omega*alpha*cos(omega*current_t+i*delta);
       ddref_phi(i)=-omega^2*alpha*sin(omega*current_t+i*delta);
    end
    
%计算控制信号
%    u_bar=ddref_phi+K_D*(dref_phi-last_dphi)+K_P*(ref_phi-last_phi);%运动控制（无轨迹跟踪）
      u_bar=K_P*(ref_phi-last_phi)-K_D*last_dphi; %轨迹跟踪控制器
     
%求解蛇形机器人模型微分方程  
    [t3, out]=ode45('ode_snakeRob',interval,init_config_vars);
    
%更新状态向量
    init_config_vars=out(end,:);

 % phi1：phi9
    temp_phi=out(end,1:1:9);
    phi1(k,1)=temp_phi(1); phi2(k,1)=temp_phi(2);phi3(k,1)=temp_phi(3); phi4(k,1)=temp_phi(4); phi5(k,1)=temp_phi(5);
    phi6(k,1)=temp_phi(6); phi7(k,1)=temp_phi(7); phi8(k,1)=temp_phi(8); phi9(k,1)=temp_phi(9);
 % Delta_theta
    Delta_theta(k,1)=delta_theta;
 % theta10
    theta10(k,1)=out(end,10);
 %px,py
    temp_p=out(end,11:1:12);
    px(k,1)=temp_p(1); py(k,1)=temp_p(2);
 %phidt1：phidt10
    temp_d_phi=out(end,13:1:21);
    d_phi1(k,1)=temp_d_phi(1); d_phi2(k,1)=temp_d_phi(2);d_phi3(k,1)=temp_d_phi(3); d_phi4(k,1)=temp_d_phi(4); d_phi5(k,1)=temp_d_phi(5);
    d_phi6(k,1)=temp_d_phi(6); d_phi7(k,1)=temp_d_phi(7); d_phi8(k,1)=temp_d_phi(8); d_phi9(k,1)=temp_d_phi(9);
 %theta10dt
    d_theta10(k,1)=out(end,10);
 %pxdt,pydt
    temp_d_p=out(end,11:1:12);
    d_px(k,1)=temp_d_p(1); d_py(k,1)=temp_d_p(2);
    Time(k,1)=current_t;
    ref_phi1(k,1)=ref_phi(1);
 %P_x, P_y
     P_x(k,1)=last_px;
     P_y(k,1)=last_py;
 %Phi0
     Phi0(k,1)=phi0;
     Phi1(k,1)=init_config_vars(1);
     Phi2(k,1)=init_config_vars(2);
     Phi3(k,1)=init_config_vars(3);
     V_x(k,1)=init_config_vars(23);
     V_y(k,1)=init_config_vars(24);
     V(k,1)=((init_config_vars(23))^2+(init_config_vars(24))^2)^0.5;
     
      global Torque
     Energy(k,1)=(Torque'*Torque);
     Velocity(k,1)=((d_px(k,1))^2+(d_py(k,1)^2))^0.5;
     Efficiency(k,1)=Velocity(k,1)/Energy(k,1);
%      if Efficiency(k,1)>30
%          Efficiency(k,1)=30-rand;
%      end
     C(k,1)=sum(Efficiency(1:k));
     sum_py=sum_py+abs(last_py);
     Sum_py(k,1)=sum_py;
end

% 画图
time=linspace(t0,tfin,num_seg_t);
  %质心位置
% figure 
% subplot(2,1,1)
% plot(P_x,P_y)
% legend('CM position')
% grid on
% subplot(2,1,2);
% plot(Time,V)
% legend('V')
% grid on

% figure 
% subplot(4,1,1)
% plot(P_x,P_y)
% legend('Position of CM')
% grid on
% subplot(4,1,2)
% plot(Time,P_x)
% legend('P_x')
% grid on
% subplot(4,1,3)
% plot(Time,P_y)
% legend('P_y')
% grid on
% subplot(4,1,4)
% plot(Time,V)
% legend('velocity')
% grid on

figure
plot(P_x,P_y,P_x,zeros(600,1))
legend('质心位置')
grid on

figure
plot(Time,Sum_py)
legend('偏差累计和')
grid on

figure
plot(Time,Efficiency)
legend('效率')
grid on

figure
plot(Time,C)
legend('运动效率之和')
grid on
%%