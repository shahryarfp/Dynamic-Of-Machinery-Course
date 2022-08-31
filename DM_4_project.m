clc;clear;
syms t4 t9 t6 t3 t10;
R7=10;R8=15;R4=48.270045;R3=49.329413;R10=10;
R9=30;R6=20;O2O6=55;
DFE=pi/2;CO2B=pi/2;
t7=5*pi/4;t8=t7+pi/2; %t = theta
w2=1;
a2=0; %a+number => a means alpha
VFF=[];aFF=[];VDD=[];aDD=[];VEE=[];aEE=[]; %for plot

T2=[];FF06=[];FF56=[];FF3D=[];FF4E=[];FFB3=[];FFC4=[];FFO2=[];
for time= linspace(0,2*pi,200)
    t7=w2*time+5*pi/4;t8=t7+pi/2;
    
    %position
    F=@(x)[R7*cos(t7)+R4*cos(x(1))+R9*cos(x(2))+R6*cos(x(3))+O2O6*cos(3*pi/2)
        R8*cos(t8)+R3*cos(x(4))+R10*cos(x(5))+R6*cos(x(3))+O2O6*cos(3*pi/2)
        R7*sin(t7)+R4*sin(x(1))+R9*sin(x(2))+R6*sin(x(3))+O2O6*sin(3*pi/2)
        R8*sin(t8)+R3*sin(x(4))+R10*sin(x(5))+R6*sin(x(3))+O2O6*sin(3*pi/2)
        x(2)-x(5)-pi/2];
    x0=[70,150,45,130,80]*pi/180;
    x=fsolve(F,x0);
    t4=x(1);t9=x(2);t6=x(3);t3=x(4);t10=x(5);
    
    %velocity
    G=@(x)[-R7*w2*sin(t7)-R4*x(1)*sin(t4)-R9*x(2)*sin(t9)-R6*x(3)*sin(t6)
      -R8*w2*sin(t8)-R3*x(4)*sin(t3)-R10*x(5)*sin(t10)-R6*x(3)*sin(t6)
      R7*w2*cos(t7)+R4*x(1)*cos(t4)+R9*x(2)*cos(t9)+R6*x(3)*cos(t6)
      R8*w2*cos(t8)+R3*x(4)*cos(t3)+R10*x(5)*cos(t10)+R6*x(3)*cos(t6)
      x(2)-x(5)];
    x0=[.5,.5,.5,.5,.5];
    x=fsolve(G,x0);
    w4=x(1);w9=x(2);w6=x(3);w3=x(4);w10=x(5);
    
    %acceleration
    H=@(x)[R7*a2*sin(t7)+R7*w2^2*cos(t7)+ R4*x(1)*sin(t4)+R4*w4^2*cos(t4)+ R9*x(2)*sin(t9)+R9*w9^2*cos(t9)+ R6*x(3)*sin(t6)+R6*w6^2*cos(t6)
      R8*a2*sin(t8)+R8*w2^2*cos(t8)+ R3*x(4)*sin(t3)+R3*w3^2*cos(t3)+ R10*x(5)*sin(t10)+R10*w10^2*cos(t10)+ R6*x(3)*sin(t6)+ R6*w6^2*cos(t6)
      R7*a2*cos(t7)-R7*w2^2*sin(t7)+R4*x(1)*cos(t4)-R4*w4^2*sin(t4)+R9*x(2)*cos(t9)-R9*w9^2*sin(t9)+R6*x(3)*cos(t6)-R6*w6^2*sin(t6)
      R8*a2*cos(t8)-R8*w2^2*sin(t8)+R3*x(4)*cos(t3)-R3*w3^2*sin(t3)+R10*x(5)*cos(t10)-R10*w10^2*sin(t10)+R6*x(3)*cos(t6)-R6*w6^2*sin(t6)
      x(2)-x(5)];
    x0=[.5,.5,.5,.5,.5];
    x=fsolve(H,x0);
    a4=x(1);a9=x(2);a6=x(3);a3=x(4);a10=x(5);
    
    VF=R6*w6*[cos(t6+3*pi/2), sin(t6+3*pi/2)];
    aF=R6*a6*[cos(t6+3*pi/2), sin(t6+3*pi/2)]+R6*w6^2*[cos(t6), sin(t6)];
    VD=VF+R10*w10*[cos(t10+3*pi/2), sin(t10+3*pi/2)];
    aD=aF+R10*w10^2*[cos(t10), sin(t10)]+R10*a10*[cos(t10+3*pi/2), sin(t10+3*pi/2)];
    VE=VF+R9*w9*[cos(t9+3*pi/2), sin(t9+3*pi/2)];
    aE=aF+R9*w9^2*[cos(t9), sin(t9)]+R9*a9*[cos(t9+3*pi/2), sin(t9+3*pi/2)];
    aG6=R6/2*a6*[cos(t6+3*pi/2), sin(t6+3*pi/2)] + R6/2*w6^2*[cos(t6), sin(t6)];
    C=R7*cos(t7);C=[C R7*sin(t7)];
    B=R8*cos(t8);B=[B R8*sin(t8)];
    E=C(1)+R4*cos(t4);E=[E C(2)+R4*sin(t4)];
    D=B(1)+R3*cos(t3);D=[D B(2)+R3*sin(t3)];
    F=E(1)+R9*cos(t9);F=[F E(2)+R9*sin(t9)];
    RG5=(F+E+D)/3-F;
    tG5=atan(3);
    aG5=aF+norm(RG5)*a10*[cos(t10+3*pi/2+tG5), sin(t10+3*pi/2+tG5)]+norm(RG5)*w10^2*[cos(t10+tG5), sin(t10+tG5)];
    aG3=aD+R3/2*w3^2*[cos(t3), sin(t3)]+R3/2*a3*[cos(t3+3*pi/2), sin(t3+3*pi/2)];
    aG4=aE+R4/2*w4^2*[cos(t4), sin(t4)]+R4/2*a4*[cos(t4+3*pi/2), sin(t4+3*pi/2)];
    RG2=(B+C)/3;
    tG2=atan(1.5);
    aG2=norm(RG2)*w2^2*[cos(t7+tG2), sin(t7+tG2)];
    
    R60=[0 55]-([0 55]+F)/2;R60=R60/100;
    R6F=F-([0 55]+F)/2;R6F=R6F/100;
    R5F=-RG5;R5F=R5F/100;
    R5D=D-(F+E+D)/3;R5D=R5D/100;
    R5E=E-(F+E+D)/3;R5E=R5E/100;
    R3D=D-(B+D)/2;R3D=R3D/100;
    R3B=B-(B+D)/2;R3B=R3B/100;
    R4E=E-(E+C)/2;R4E=R4E/100;
    R4C=C-(E+C)/2;R4C=R4C/100;
    R2O=-(B+C)/3;R2O=R2O/100;
    R2C=C-(B+C)/3;R2C=R2C/100;
    R2B=B-(B+C)/3;R2B=R2B/100;
    T6=20; %N.m
    
    m6=0.2;
    m5=0.015;
    m2=0.0075;
    m3=0.4933;
    m4=0.4827;
    
    I3=0.01;
    I6=0.00067;
    I4=0.00937;
    I2=4.0625e-5;
    I5=205e-4;
    a5=a9;
    
    l1=[1,0,1,0,0,0,0,0,0,0,0,0,0,0,0];
    l2=[0,1,0,1,0,0,0,0,0,0,0,0,0,0,0];
    l3=[-R60(2),R60(1),-R6F(2),R6F(1),0,0,0,0,0,0,0,0,0,0,0];
    l4=[0,0,-1,0,1,0,1,0,0,0,0,0,0,0,0];
    l5=[0,0,0,-1,0,1,0,1,0,0,0,0,0,0,0];
    l6=[0,0,R5F(2),-R5F(1),-R5D(2),R5D(1),-R5E(2),R5E(1),0,0,0,0,0,0,0];
    l7=[0,0,0,0,-1,0,0,0,1,0,0,0,0,0,0];
    l8=[0,0,0,0,0,-1,0,0,0,1,0,0,0,0,0];
    l9=[0,0,0,0,R3D(2),-R3D(1),0,0,1,-R3B(2),R3B(1),0,0,0,0];
    l10=[0,0,0,0,0,0,-1,0,0,0,1,0,0,0,0];
    l11=[0,0,0,0,0,0,0,-1,0,0,0,1,0,0,0];
    l12=[0,0,0,0,0,0,R4E(2),-R4E(1),0,0,-R4C(2),R4C(1),0,0,0];
    l13=[0,0,0,0,0,0,0,0,-1,0,-1,0,1,0,0];
    l14=[0,0,0,0,0,0,0,0,0,-1,0,-1,0,1,0];
    l15=[0,0,0,0,0,0,0,0,R2B(2),-R2B(1),R2C(2),-R2C(1),-R2O(2),R2O(1),1];
    
    known=[m6*aG6(1)/100;m6*aG6(2)/100;I6*a6-T6;
        m5*aG5(1)/100;m5*aG5(2)/100-100;I5*a5;
        m3*aG3(1)/100-100;m3*aG3(2)/100;I3*a3;
        m4*aG4(1)/100-100;m4*aG4(2)/100;I4*a4;
        m2*aG2(1)/100;m2*aG2(2)/100;I2*a2];
    
    Ans=inv([l1;l2;l3;l4;l5;l6;l7;l8;l9;l10;l11;l12;l13;l14;l15])*known;
    F06=Ans(1);F06=[F06,Ans(2)];FF06=[FF06 norm(F06)];
    F56=Ans(3);F56=[F56,Ans(4)];FF56=[FF56 norm(F56)];
    F3D=Ans(5);F3D=[F3D,Ans(6)];FF3D=[FF3D norm(F3D)];
    F4E=Ans(7);F4E=[F4E,Ans(8)];FF4E=[FF4E norm(F4E)];
    FB3=Ans(9);FB3=[FB3,Ans(10)];FFB3=[FFB3 norm(FB3)];
    FC4=Ans(11);FC4=[FC4,Ans(12)];FFC4=[FFC4 norm(FC4)];
    FO2=Ans(13);FO2=[FO2,Ans(14)];FFO2=[FFO2 norm(FO2)];
    T2=[T2 Ans(15)];
end

%plot T2
subplot(1,2,1);
Time= linspace(0,2*pi,200);
plot(Time,T2,'b');
grid on;
title('T2-time');
ylim([-20,20]);
xlim([0,2*pi]);
ylabel('T2(N.m)');
xlabel('time(s)');

%plot forces
subplot(1,2,2);
plot(Time,FF06,'r');
grid on;
title('F-time');
ylim([0,300]);
xlim([0,2*pi]);
ylabel('F(N)');
xlabel('time(s)');
hold on;
plot(Time,FF56,'b');
plot(Time,FF3D,':k');
plot(Time,FF4E,'--');
plot(Time,FFB3,'*');
plot(Time,FFC4,'b');
plot(Time,FFO2,'.');
legend('F-O6','F-F','F-D','F-E','F-B','F-C','F-O2');