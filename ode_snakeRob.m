function f=ode_snakeRob(t,inp)
% phi1=inp(1);phi2=inp(2);phi3=inp(3);phi4=inp(4);phi5=inp(5);phi6=inp(6);phi7=inp(7);phi8=inp(8);phi9=inp(9);
% theta10=inp(10);px=inp(11);py=inp(12);
% dphi1=inp(13);dphi2=inp(14);dphi3=inp(15);dphi4=inp(16);dphi5=inp(17);dphi6=inp(18);dphi7=inp(19);dphi8=inp(20);dphi9=inp(21);
% dtheta10=inp(22);dpx=inp(23);dpy=inp(24);
global N m l j cn ct 
x(1)=inp(1);x(2)=inp(2);x(3)=inp(3);x(4)=inp(4);x(5)=inp(5);x(6)=inp(6);x(7)=inp(7);x(8)=inp(8);x(9)=inp(9);x(10)=inp(10);x(11)=inp(11);x(12)=inp(12);
x(13)=inp(13);x(14)=inp(14);x(15)=inp(15);x(16)=inp(16);x(17)=inp(17);x(18)=inp(18);x(19)=inp(19);x(20)=inp(20);x(21)=inp(21);x(22)=inp(22);x(23)=inp(23);x(24)=inp(24);

N=10;m=1;l=0.14;ct=1;cn=3;j=0.0016;
A=[1,1,0,0,0,0,0,0,0,0
   0,1,1,0,0,0,0,0,0,0
   0,0,1,1,0,0,0,0,0,0
   0,0,0,1,1,0,0,0,0,-1
   0,0,0,0,1,1,0,0,0,0
   0,0,0,0,0,1,1,0,0,0
   0,0,0,0,0,0,1,1,0,0
   0,0,0,0,0,0,0,1,1,0
   0,0,0,0,0,0,0,0,1,1];%% size(A)=9*10
D=[1,-1,0,0,0,0,0,0,0,0
   0,1,-1,0,0,0,0,0,0,0
   0,0,1,-1,0,0,0,0,0,0
   0,0,0,1,-1,0,0,0,0,0
   0,0,0,0,1,-1,0,0,0,0
   0,0,0,0,0,1,-1,0,0,0
   0,0,0,0,0,0,1,-1,0,0
   0,0,0,0,0,0,0,1,-1,0
   0,0,0,0,0,0,0,0,1,-1];%% size(D)=9*10
H=[1,1,1,1,1,1,1,1,1,1
   0,1,1,1,1,1,1,1,1,1
   0,0,1,1,1,1,1,1,1,1
   0,0,0,1,1,1,1,1,1,1
   0,0,0,0,1,1,1,1,1,1
   0,0,0,0,0,1,1,1,1,1
   0,0,0,0,0,0,1,1,1,1
   0,0,0,0,0,0,0,1,1,1
   0,0,0,0,0,0,0,0,1,1
   0,0,0,0,0,0,0,0,0,1];%% size(H)=10*10
E1=ones(N,1); %% size(E1)=10*1
E2=[E1,zeros(N,1);zeros(N,1),E1]; %%size(E1)=20*2
In=eye(N); %% size(In)=10*10
% Theta=zeros(N,1);
Theta=H*[x(1);x(2);x(3);x(4);x(5);x(6);x(7);x(8);x(9);x(10)];%% size(Theta)=10*1
% Thetadt=zeros(N,1);
Thetadt=H*[x(13);x(14);x(15);x(16);x(17);x(18);x(19);x(20);x(21);x(22)];%% size(Thetadt)=10*1
% theta1=0;theta2=0;theta3=0;theta4=0;theta5=0;theta6=0;theta7=0;theta8=0;theta9=0;theta10=0;
theta1=Theta(1,1);theta2=Theta(2,1);theta3=Theta(3,1);theta4=Theta(4,1);theta5=Theta(5,1);theta6=Theta(6,1);theta7=Theta(7,1);theta8=Theta(8,1);theta9=Theta(9,1);theta10=Theta(10,1);
% Sintheta=zeros(10,1);
% Costheta=zeros(10,1);
Sintheta=[sin(theta1);sin(theta2);sin(theta3);sin(theta4);sin(theta5);sin(theta6);sin(theta7);sin(theta8);sin(theta9);sin(theta10)];
Costheta=[cos(theta1);cos(theta2);cos(theta3);cos(theta4);cos(theta5);cos(theta6);cos(theta7);cos(theta8);cos(theta9);cos(theta10)];
%%size(Sintheta)=10*1
%%size(Costheta)=10*1
% Stheta=zeros(10,10);
% Ctheta=zeros(10,10);
Stheta=diag(Sintheta);%%size(Stheta)=10*10
Ctheta=diag(Costheta);%%size(Ctheta)=10*10
V=(A'/(D*D'))*A;%%size(V)=10*10
K=(A'/(D*D'))*D;%%size(K)=10*10
% Mtheta=zeros(10,10);
% Mbar=zeros(12,12);
Mtheta=j*In+m*l^2*Stheta*V*Stheta+m*l^2*Ctheta*V*Ctheta;%%size(Mtheta)=10*10
Mbar=[H'*Mtheta*H,zeros(N,2);zeros(2,N),N*m*eye(2)];%% size(Mbar)=12*12
M11bar=Mbar(1:N-1,1:N-1);%%size(M11)=(N-1)*(N-1)
M12bar=Mbar(1:N-1,N:N+2);%%size(M12)=(N-1)*(3)
M21bar=Mbar(N:N+2,1:N-1);%%size(W21)=3*9
M22bar=Mbar(N:N+2,N:N+2);%%size(M22)=3*3
Phibardt=zeros(10,1);
Wphibar=zeros(10,1);
% Wbar=zeros(12,1);
%%Gbar=zeros(12,20)
Phibardt=[x(13);x(14);x(15);x(16);x(17);x(18);x(19);x(20);x(21);x(22)];%%size(Phibar)=10*1
Wphibar=m*l^2*Stheta*V*Stheta-m*l^2*Ctheta*V*Ctheta;%%size(Wphibar)=10*1
Wbar=[H'*Wphibar*diag(H*Phibardt)*H*Phibardt;zeros(2,1)];%%size(Wbar)=12*1
Gbar=[-l*H'*Stheta*K,l*H'*Ctheta*K;-E1',zeros(1,N);zeros(1,N),-E1'];%%size(Gbar)=12*20
W1bar=Wbar(1:9,1);%%size(W1bar£©=9*1
W2bar=Wbar(10:12,1);%%size(W2bar)=3*1
G1bar=Gbar(1:9,:); %% size(G1bar)=9*20
G2bar=Gbar(10:12,:);%%size(G2bar)=3*20
% Xdt=zeros(10,1);
% Ydt=zeros(10,1);
Fr=zeros(20,1);
Xdt=l*K'*Stheta*Thetadt+E1*x(23);%%size(Xdt)=10*1
Ydt=-l*K'*Ctheta*Thetadt+E1*x(24);%%size(Ydt)=10*1
Fr=-[ct*(Ctheta)^2+cn*(Stheta)^2,(ct-cn)*Stheta*Ctheta;(ct-cn)*Stheta*Ctheta,ct*(Stheta)^2+cn*(Ctheta)^2]*[Xdt;Ydt];%%size(Fr)=20*1
%%A1=-(W2bar+G2bar*Fr)/M22bar;%%size(A)=3*1
%%B=-M21bar/M22bar;%%size(B)=3*9
% C=zeros(3,1);
% C1=zeros(3,9);
% C2=zeros(3,1);
A8=G2bar*Fr;
% A9=A8/M22bar;
C=-M22bar^-1*(W2bar+G2bar*Fr);%%size(C)=3*1
C1=-((M22bar)^-1)*M21bar;%%size(C1)=3*9

global u_bar

C2=C+C1*[u_bar(1);u_bar(2);u_bar(3);u_bar(4);u_bar(5);u_bar(6);u_bar(7);u_bar(8);u_bar(9)];%%size(C2)=3*1
dx1=x(13);dx2=x(14);dx3=x(15);dx4=x(16);dx5=x(17);dx6=x(18);dx7=x(19);dx8=x(20);dx9=x(21);
dx10 =x(22);
dx11 =x(23);
dx12 =x(24);
dx13=u_bar(1,1);dx14=u_bar(2,1);dx15=u_bar(3,1);dx16=u_bar(4,1);dx17=u_bar(5,1);dx18=u_bar(6,1);dx19=u_bar(7,1);dx20=u_bar(8,1);dx21=u_bar(9,1);
dx22=C2(1,1);
dx23=C2(2,1);
dx24=C2(3,1);

global Torque
Torque=(M11bar-(M12bar/(M22bar))*M21bar)*u_bar+W1bar+G1bar*Fr-(M12bar/(M22bar))*(W2bar+G2bar*Fr);
f =[dx1;dx2;dx3;dx4;dx5;dx6;dx7;dx8;dx9;dx10;dx11;dx12;dx13;dx14;dx15;dx16;dx17;dx18;dx19;dx20;dx21;dx22;dx23;dx24];
% f1=[dphi1;dphi2;dphi3;dphi4;dphi5;dphi6;dphi7;dphi8;dphi9;dtheta10;dpx;dpy;ddphi1;ddphi2;ddphi3;ddphi4;ddphi5;ddphi6;ddphi7;ddphi8;ddphi9;ddtheta10;ddpx;ddpy]
