a=100;
s=2;
u=10;
xd=1:20;
yd=a*exp(-xd*s);
func = @(x,xdata)x(1)*exp(x(2)*(xdata-x(3)));
x0=[2,30];
P = lsqcurvefit(func,x0,xd,yd)