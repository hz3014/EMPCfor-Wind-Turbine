function x=precompute(x1,u)
global R
global Ts
global T
global Vt
global M

J=1/9*R^2*M;
B=0.01;
density=1.225;
zeta=0.131;
x(1)=x1;

for i=1:T-1
    lembda=x(i)*R/Vt;  %>2.5, 5
    lembdai=1/(1/(lembda+0.08*zeta)-0.035/(zeta^3+1)); 
    Cp=0.22*(116/lembdai-0.4*zeta-5)*exp(-12.5/lembdai);
    Ta=density*pi*R^2*Vt^3*Cp/(2*x(i));
    %W.t(i)=Ta/J-u(i)/J;
    W.t(i)=Ta/J-B*x(i)/J-u(i)/J;
    x(i+1)=x(i)+Ts*W.t(i);
end;