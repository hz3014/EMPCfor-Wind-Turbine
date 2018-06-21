function [c,ceq]= nonlcon(x)
global R
global Ts
global T
global Vt
global M

J=1/9*R^2*M;
B=0.01;
density=1.225;
zeta=0.131;
ceq=[];
c=[];

for i=1:T-1
    lembda=x(i)*R/Vt;  %>2.5, 5
    lembdai=1/(1/(lembda+0.08*zeta)-0.035/(zeta^3+1)); 
    Cp=0.22*(116/lembdai-0.4*zeta-5)*exp(-12.5/lembdai);
    Ta=density*pi*R^2*Vt^3*Cp/(2*x(i));
    %W.t(i)=Ta/J-x(i+T)/J;
    W.t(i)=Ta/J-B*x(i)/J-x(i+T)/J;
    %W.t(i)=Ta/J-B*x(i)/J-6500/J;
    %x(i+1)=x(i)+Ts*W.t
    ceq=[ceq;x(i+1)-Ts*W.t(i)-x(i)];
    %ceq=[W.t(T-1)]
end
 
end