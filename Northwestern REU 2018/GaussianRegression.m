a=100;
s=2;
u=10;
x=1:20;
hold on
y=a*exp(-(1/2)*((x-u)/s).^2)./(s*sqrt(2*pi));
plot(x,y,'ro')
y1=log(y);
P=polyfit(x,y1,2)
std=1/sqrt(-2*P(1))
mean=P(2)*std*std
A=max(y)*sqrt(2*pi)*std
plot(x,(A/(std*sqrt(2*pi)))*exp((-1/2)*((x-mean)/std).^(2)),'k')
hold off