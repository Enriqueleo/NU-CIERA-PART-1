import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import glob
import errno
#I read in the file of the Type Ia Supernova
#Fitting model equation
def f(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))
F=np.loadtxt\
('C:/Users/enriq/Desktop/Northwestern REU 2018/DATASET/sub DATA 1/SN2007asopt_gFilter_.txt')
#The three col of information from the data
t=F[:,0]
mag=F[:,1]
err=F[:,2]
#Finding the parameters of the model
t0=t[np.where(mag==max(mag))[0]]
bopt=(mag[len(t)-1]-mag[len(t)-3])/(t[len(t)-1]-t[len(t)-3])/max(mag)
coef,covar=spo.curve_fit(f,t,mag,[max(mag),-bopt,-0.01,0.01,0.28,20,t0],sigma=err)
pa,covar=spo.curve_fit(f,t,mag,[coef[0],coef[1],coef[2],coef[3],coef[4],coef[5],coef[6]],sigma=err)
x=np.linspace(min(t),max(t)+0.5,81)
plt.figure(figsize=(6,6))
plt.plot(x,f(x,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),color='purple')
plt.errorbar(t,mag,yerr=err,fmt='o',color='red')
plt.ylabel('Magnitude',Fontsize='16',color='black')
plt.xlabel('Epoch (days)',Fontsize='16',color='Black')
plt.title("Light Curve",Fontsize='18',color='Black')

chisq=sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7)
print(f'𝛘2={chisq:5.4}')
p=['A','b','c','d','k','σ','t0']
for i in range(len(pa)):
    print(f'{p[i]}={pa[i]:4.6}±{np.diag(covar)[i]:4.3}')