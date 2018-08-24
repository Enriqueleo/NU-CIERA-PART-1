import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd

def constant(A,c,d):
        return (A-A*d)/(1+c)
def exp1(x,c,b,s):
        return 1-b*x+c*np.exp(-(x/s)**2)
def exp2(x,d,k):
        return 1-d*np.exp(-k*x)

def fitting(x,A,b,c,d,k,s,t):
        return (constant(A,c,d)*exp1(x-t,c,b,s))/exp2(x-t,d,k)

x=np.array(range(int(53773),int(53900)+1))
y=fitting(x,-14.802,(-852)/(10**6),-0.116,0.00123,0.289,21.4,53788.48)
plt.plot(x,y)
print(fitting(53788.48,-14.802,(-852)/(10**6),-0.116,0.00123,0.289,21.4,53788.48))