clc;clear;
syms t4 t9 t6 t3 t10;
R7=10;R8=15;R4=48.270045;R3=49.329413;R10=10;
R9=30;R6=20;O2O6=55;
DFE=pi/2;CO2B=pi/2;
t7=5*pi/4;t8=t7+pi/2; %t = theta
w2=1;
a2=0; %a+number => a means alpha
VFF=[];aFF=[];VDD=[];aDD=[];VEE=[];aEE=[]; %for plot

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
    
    G=@(x)[-R7*w2*sin(t7)-R4*x(1)*sin(t4)-R9*x(2)*sin(t9)-R6*x(3)*sin(t6)
      -R8*w2*sin(t8)-R3*x(4)*sin(t3)-R10*x(5)*sin(t10)-R6*x(3)*sin(t6)
      R7*w2*cos(t7)+R4*x(1)*cos(t4)+R9*x(2)*cos(t9)+R6*x(3)*cos(t6)
      R8*w2*cos(t8)+R3*x(4)*cos(t3)+R10*x(5)*cos(t10)+R6*x(3)*cos(t6)
      x(2)-x(5)];
    x0=[.5,.5,.5,.5,.5];
    x=fsolve(G,x0);
    w4=x(1);w9=x(2);w6=x(3);w3=x(4);w10=x(5);
    
    %velocity
    H=@(x)[R7*a2*sin(t7)+R7*w2^2*cos(t7)+ R4*x(1)*sin(t4)+R4*w4^2*cos(t4)+ R9*x(2)*sin(t9)+R9*w9^2*cos(t9)+ R6*x(3)*sin(t6)+R6*w6^2*cos(t6)
      R8*a2*sin(t8)+R8*w2^2*cos(t8)+ R3*x(4)*sin(t3)+R3*w3^2*cos(t3)+ R10*x(5)*sin(t10)+R10*w10^2*cos(t10)+ R6*x(3)*sin(t6)+ R6*w6^2*cos(t6)
      R7*a2*cos(t7)-R7*w2^2*sin(t7)+R4*x(1)*cos(t4)-R4*w4^2*sin(t4)+R9*x(2)*cos(t9)-R9*w9^2*sin(t9)+R6*x(3)*cos(t6)-R6*w6^2*sin(t6)
      R8*a2*cos(t8)-R8*w2^2*sin(t8)+R3*x(4)*cos(t3)-R3*w3^2*sin(t3)+R10*x(5)*cos(t10)-R10*w10^2*sin(t10)+R6*x(3)*cos(t6)-R6*w6^2*sin(t6)
      x(2)-x(5)];
    x0=[.5,.5,.5,.5,.5];
    x=fsolve(H,x0);
    a4=x(1);a9=x(2);a6=x(3);a3=x(4);a10=x(5);
    
    %acceleration
    VF=R6*w6*[cos(t6+3*pi/2), sin(t6+3*pi/2)];
    aF=R6*a6*[cos(t6+3*pi/2), sin(t6+3*pi/2)]+R6*w6^2*[cos(t6), sin(t6)];
    VD=VF+R10*w10*[cos(t10+3*pi/2), sin(t10+3*pi/2)];
    aD=aF+R10*w10^2*[cos(t10), sin(t10)]+R10*a10*[cos(t10+3*pi/2), sin(t10+3*pi/2)];
    VE=VF+R9*w9*[cos(t9+3*pi/2), sin(t9+3*pi/2)];
    aE=aF+R9*w9^2*[cos(t9), sin(t9)]+R9*a9*[cos(t9+3*pi/2), sin(t9+3*pi/2)];
    VFF=[VFF norm(VF)];
    aFF=[aFF norm(aF)];
    VDD=[VDD norm(VD)];
    aDD=[aDD norm(aD)];
    VEE=[VEE norm(VE)];
    aEE=[aEE norm(aE)];
end

%plot velocity
subplot(1,2,1);
Time= linspace(0,2*pi,200);
plot(Time,VFF,'b');
grid on;
title('V-time');
ylim([-2,20]);
xlim([0,2*pi]);
ylabel('V(cm/s)');
xlabel('time(s)');
hold on;
plot(Time,VDD,'r');
plot(Time,VEE,'k');
legend('VF','VD','VE');

%plot acceleration
subplot(1,2,2);
plot(Time,aFF,'b');
grid on;
title('A-time');
ylim([-2,40]);
xlim([0,2*pi]);
ylabel('A(cm/s^2)');
xlabel('time(s)');
hold on;
plot(Time,aDD,'r');
plot(Time,aEE,'k');
legend('AF','AD','AE');