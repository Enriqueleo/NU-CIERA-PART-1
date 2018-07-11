x=-1:6;
y=2*x+10;
Error=ones(size(y));
%
delta=sum(Error.^(-2).*x.^2)*sum((Error.^(-2)))-sum(Error.^(-2).*x).^2;
Sum1=sum(Error.^(-2).*y);
Sum2=sum(Error.^(-2).*x);
Sum3=sum(Error.^(-2).*x.*y);
slope=(sum((Error.^(-2)))*Sum3-Sum1*Sum2)/delta;
yinter=(sum(Error.^(-2).*x.^2)*Sum1-Sum2*Sum3)/delta;
P=polyfit(x,y,1);
hold on
plot(x,P(1)*x+P(2),'r')
errorbar(x,y,Error,'ko');
hold off
fprintf('%4.1f is the slope and %4.1f y-intercept\n',P(1),P(2))
