clc;clear;
syms t4 t9 t6 t3 t10;
R7=10;R8=15;R4=48.270045;R3=49.329413;R10=10;
R9=30;R6=20;O2O6=55;
DFE=pi/2;CO2B=pi/2;
t7=5*pi/4;t8=t7+pi/2; %t = theta
w2=1;
a2=0; %in a+number => a means alpha

time_step=0.1;
for i=1:360
    t7=t7+time_step*w2;t8=t7+pi/2;
    G=@(x)[R7*cos(t7)+R4*cos(x(1))+R9*cos(x(2))+R6*cos(x(3))+O2O6*cos(3*pi/2)
        R8*cos(t8)+R3*cos(x(4))+R10*cos(x(5))+R6*cos(x(3))+O2O6*cos(3*pi/2)
        R7*sin(t7)+R4*sin(x(1))+R9*sin(x(2))+R6*sin(x(3))+O2O6*sin(3*pi/2)
        R8*sin(t8)+R3*sin(x(4))+R10*sin(x(5))+R6*sin(x(3))+O2O6*sin(3*pi/2)
        x(2)-x(5)-pi/2];
    x0=[70,150,45,130,80]*pi/180;
    x=fsolve(G,x0);
    t4=x(1);t9=x(2);t6=x(3);t3=x(4);t10=x(5);
    C=R7*cos(t7);C=[C R7*sin(t7)];
    B=R8*cos(t8);B=[B R8*sin(t8)];
    E=C(1)+R4*cos(t4);E=[E C(2)+R4*sin(t4)];
    D=B(1)+R3*cos(t3);D=[D B(2)+R3*sin(t3)];
    F=E(1)+R9*cos(t9);F=[F E(2)+R9*sin(t9)];
    
    plot([0;C(1)],[0;C(2)],'b');
    hold on;
    xlim([-30 30]);
    ylim([-20 75]);
    plot([0;B(1)],[0;B(2)],'b');
    plot([0;-3],[0;-3],'b');
    plot([0;3],[0;-3],'b');
    plot([3;-3],[-3;-3],'b');
    
    plot([0;-3],[0+55;-3+55],'b');
    plot([0;3],[0+55;-3+55],'b');
    plot([3;-3],[-3+55;-3+55],'b');
    plot([B(1);C(1)],[B(2);C(2)],'b');
    plot([B(1);C(1)],[B(2);C(2)],'b');
    plot([B(1);D(1)],[B(2);D(2)],'b');
    plot([C(1);E(1)],[C(2);E(2)],'b');
    plot([0;F(1)],[55;F(2)],'b');
    plot([F(1);E(1)],[F(2);E(2)],'b');
    plot([F(1);D(1)],[F(2);D(2)],'b');
    plot([D(1);E(1)],[D(2);E(2)],'b');
    
    pause(0.001);
    hold off;
end