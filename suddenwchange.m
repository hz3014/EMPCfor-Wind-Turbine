clear;
global V
global T
global Vt
global Ts
global R
global M
R=76;
Vt=9;%initial wind speed setting
T=6; 
Ts=0.1;
M=4800;
B=0.01;
IP=0;
density=1.225;
zeta=0.131;
J=1/9*R^2*M;
lb1=[0.1;0];
ub1=[0.45;8000000];
y0=[0.2;2000000]; 
A1=[];
b1=[];
Aeq1=[];
beq1=[];

[y,fvall,exitflag1,output1]=fmincon('powerextracted',y0,A1,b1,Aeq1,beq1,lb1,ub1,'dude')
Xs=y(1);
X0=y(2);
x1=0.45;%initial rotor speed setting
x(1)=x1;
u(1:T)=7000000;
y=precompute(x1,u);
lb(1)=x(1);lb(2:T)=0.1;lb(T+1:2*T-1)=0;
ub(1)=x(1);ub(2:T)=0.45;ub(T+1:2*T-1)=8000000 ;
x0(1)=x(1);x0(2:3)=x(1);x0(4:5)=(x(1)+Xs)/2;x0(6)=Xs;x0(T+1:2*T-1)=X0;
A=[];
b=[];
Aeq(1:T-1)=0;Aeq(T)=1;Aeq(T+1:2*T-1)=0;
beq=[Xs];

options = optimoptions('fmincon','MaxFunctionEvaluations',6000,'ConstraintTolerance',5e-3);
[x,fval,exitflag,output]=fmincon('integratedpower',x0,A,b,Aeq,beq,lb,ub,'nonlcon',options);

X(1)=x(1);
for t=1:100
    V(t)=Vt;
    if t==15
        Vt=4; %sudden wind change setting
        [y,fvall,exitflag1,output1]=fmincon('powerextracted',y0,A1,b1,Aeq1,beq1,lb1,ub1,'dude');
        Xs=y(1);
        X0=y(2);
        Aeq(1:T-1)=0;Aeq(T)=1;Aeq(T+1:2*T-1)=0;
        beq=[Xs];
    end 
    IP=IP+X(t)*x(T+1);
    Torque(t)=x(T+1);
    lembda=x(1)*R/Vt;  
    lembdai=1/(1/(lembda+0.08*zeta)-0.035/(zeta^3+1));
    Cp=0.22*(116/lembdai-0.4*zeta-5)*exp(-12.5/lembdai);
    Ta=density*pi*R^2*Vt^3*Cp/(2*x(1));
    X(t+1)=X(t)+Ts*(Ta/J-B*X(t)/J-x(T+1)/J)
    x(1)=X(t+1);
    lb(1)=x(1);lb(2:T)=0.1;lb(T+1:2*T-1)=0;
    ub(1)=x(1);ub(2:T)=0.45;ub(T+1:2*T-1)=8000000;
    x0(1)=x(1);x0(2)=x(1);x0(3:4)=(x(1)+Xs)/2;x0(5:6)=Xs;x0(T+1:2*T-1)=X0;
    [x,fvalll,exitflag2,output2]=fmincon('integratedpower',x0,A,b,Aeq,beq,lb,ub,'nonlcon',options); 
end
