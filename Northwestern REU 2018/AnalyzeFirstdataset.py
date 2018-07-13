import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import time
tp=pd.read_excel('SN2006Xopt_gFilter_.xlsx',sheet_name='1')

def ef(x,c,b,d,k,s):
        return (1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x))
def fitting(x,A,b,c,d,k,s,t):
        return ((A-A*d)/(1+c))*ef(x-t,c,b,d,k,s)

t=tp['traw']
mag=tp['Magnitude']
t0=tmp[np.where(mag==max(mag))[0]]
bopt=(mag[len(t)-1]-mag[len(t)-4])/(tmp[len(t)-1]-tmp[len(t)-4])
pa,covar=spo.curve_fit(fitting,t,mag,[max(mag),-bopt,\
 -0.01,0.001,0.3,20,t0],sigma=tp['Uncertainty'])

x=np.linspace(min(t)-3,max(t)+2,200)
plt.figure(figsize=(6,6))
plt.plot(x,fitting(x,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),'r')
plt.scatter(t,mag)
plt.ylabel('Magnitude')
plt.xlabel('Epoch (days)')
plt.title("Light Curve")

parameter=['A','b','c','d','k','σ','t0']
for i in range(len(pa)):
        print(f'{parameter[i]}={pa[i]:4.6}±{np.diag(covar)[i]:4.3}')
