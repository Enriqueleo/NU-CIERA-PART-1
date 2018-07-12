import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import random as r
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

pa,covar=spo.curve_fit(fitting,tp['traw'],tp['Magnitude'],\
[-15,-8/(10**4),-0.12,0.0013,0.3,21,53790],\
sigma=tp['Uncertainty'])

plt.figure(figsize=(6,6))
plt.plot(tp['traw'],fitting(tp['traw'],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),'r')
plt.scatter(tp['traw'],tp['Magnitude'],)
plt.ylabel('Magnitude')
plt.xlabel('Epoch (days)')
print(pa)
print(np.diag(covar))