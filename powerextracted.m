function P=powerextracted(y) 
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

lembda=y(1)*R/Vt;
lembdai=1/(1/(lembda+0.08*zeta)-0.035/(zeta^3+1));
Cp=0.22*(116/lembdai-0.4*zeta-5)*exp(-12.5/lembdai);
Ta=density*pi*R^2*Vt^3*Cp/(2*y(1));
%W.t=Ta/J-y(2)/J;
W.t=Ta/J-B*y(1)/J-y(2)/J;
%Tg=Ta-J*W.t-B*y(1);
P=-y(1)*y(2);

end
