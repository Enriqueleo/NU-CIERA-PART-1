import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import random as r
import time
tp=pd.read_excel('SN2006Xopt_gFilter_.xlsx',sheet_name='1')

#plt.scatter(tp['traw'],tp['Magnitude'],)
plt.ylabel('Magnitude')
plt.xlabel('Epoch (days)')
print(tp["traw"])
def fit(x,A,b,c,d,k,s,t):
    u=x-t
    return A*((1-d)/(1+c))*((1-b*u+c*np.exp(-(u/s)**2))/(1-d*np.exp(-k*u)))
y=fit(tp["traw"],100000,1,10,2,1,20,55000)
plt.plot(tp["traw"],y,'r')