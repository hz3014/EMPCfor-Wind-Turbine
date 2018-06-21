function ip=integratedpower(x)
ip=0;
%ip=powerextracted(0.3+0.15*rand)
global T;
for i=1:T-1
    ip=ip-x(i)*x(i+T);
end
end
