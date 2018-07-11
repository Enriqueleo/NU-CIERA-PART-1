xd=-1:6;
yd=2*xd+10;
fun = @(x,xdata)x(1)*xdata+x(2);
x0=[0,5];
P = lsqcurvefit(fun,x0,xd,yd)