import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import time
tp=pd.read_excel('SN2006Xopt_gFilter_.xlsx',sheet_name='1')

def constant(A,c,d):
        return (A-A*d)/(1+c)
def exp1(x,c,b,s):
        return 1-b*x+c*np.exp(-(x/s)**2)
def exp2(x,d,k):
        return 1-d*np.exp(-k*x)

def fitting(x,A,b,c,d,k,s,t):
        return (constant(A,c,d)*exp1(x-t,c,b,s))/exp2(x-t,d,k)

tmp=tp['traw']
mag=tp['Magnitude']
t0=tmp[np.where(mag==max(mag))[0]]
pa,covar=spo.curve_fit(fitting,tmp,mag,[max(mag),0.001,\
 -0.01,0.001,0.3,20,t0],sigma=tp['Uncertainty'])

x=np.array(range(int(min(tmp)-2),int(max(tmp)+2)))
plt.figure(figsize=(6,6))
plt.plot(x,fitting(x,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),'r')
plt.scatter(tmp,mag)
plt.ylabel('Magnitude')
plt.xlabel('Epoch (days)')

parameter=['A','b','c','d','k','sigma','to']
for i in range(len(pa)):
        print(f'{parameter[i]}={pa[i]:4.6}±{np.diag(covar)[i]:4.3}')
