import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import glob
path1 = 'C:/Users/enriq/Desktop/Northwestern REU 2018/DATASET/sub DATA 1/*.txt'
path2 = 'C:/Users/enriq/Desktop/Northwestern REU 2018/DATASET/sub DATA 2/*.txt'
path3 = 'C:/Users/enriq/Desktop/Northwestern REU 2018/DATASET/sub DATA 3/*.txt'
#I read in the file of the Type Ia Supernova
#Fitting model equation
def f(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))
files1 = glob.glob(path1)
files2 = glob.glob(path2)
files3 = glob.glob(path3)

A=[]
b=[]
c=[]
d=[]
k=[]
sig=[]
tpeak=[]
chisq=[]
#Finding the parameters of each data set read in.
for name in files1:
    F=np.loadtxt(name)
    t=F[:,0]
    mag=F[:,1]
    err=F[:,2]
    t0=t[np.where(mag==max(mag))[0]]
    bopt=(mag[len(t)-1]-mag[len(t)-3])/(t[len(t)-1]-t[len(t)-3])/max(mag)
    cf,covar=spo.curve_fit(f,t,mag,[max(mag),-bopt,-0.01,0.01,0.28,20,t0],sigma=err)
    cf,covar=spo.curve_fit(f,t,mag,[max(mag),cf[1],-0.01,0.01,0.28,20,t0],sigma=err)
    cf,covar=spo.curve_fit(f,t,mag,[max(mag),cf[1],cf[2],0.01,0.28,20,t0],sigma=err)
    pa,covar=spo.curve_fit(f,t,mag,[cf[0],cf[1],cf[2],0.01,0.28,20,t0],sigma=err)
    chisq.append(sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7))
    A.append(pa[0])
    b.append(pa[1])
    c.append(pa[2])
    d.append(pa[3])
    k.append(pa[4])
    sig.append(abs(pa[5]))
    tpeak.append(pa[6])

for name2 in files2:
    F=np.loadtxt(name2)
    t=F[:,0]
    mag=F[:,1]
    err=F[:,2]
    t0=t[np.where(mag==max(mag))[0]]
    bopt=(mag[len(t)-1]-mag[len(t)-2])/(t[len(t)-1]-t[len(t)-2])/max(mag)
    pa,covar=spo.curve_fit(f,t,mag,[max(mag),-bopt,-0.01,0.011,0.27,18,t0],sigma=err)
    chisq.append(sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7))
    A.append(pa[0])
    b.append(pa[1])
    c.append(pa[2])
    d.append(pa[3])
    k.append(pa[4])
    sig.append(abs(pa[5]))
    tpeak.append(pa[6])

for name3 in files3:
    F=np.loadtxt(name3)
    t=F[:,0]
    mag=F[:,1]
    err=F[:,2]
    t0=t[np.where(mag==max(mag))[0]]
    bopt=(mag[len(t)-1]-mag[len(t)-2])/(t[len(t)-1]-t[len(t)-2])/max(mag)
    pa,covar=spo.curve_fit(f,t,mag,[max(mag),0,-0.01,0.01,0.02,20,t0],sigma=err)
    chisq.append(sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7))
    A.append(pa[0])
    b.append(pa[1])
    c.append(pa[2])
    d.append(pa[3])
    k.append(pa[4])
    sig.append(abs(pa[5]))
    tpeak.append(pa[6])

A=np.array(A)
b=np.array(b)
c=np.array(c)
d=np.array(d)
k=np.array(k)
sig=np.array(sig)
tpeak=np.array(tpeak)
#This creates a histogram of the parameters from every data set that was 
#read in.
titles=['Amplitudes','b','c','d','k','σ','t0']
param=[A,b,c,d,k,sig,tpeak]
for i in range(len(titles)):
    fig=plt.figure(figsize=(13,12))
    plt.hist(param[i],bins=np.linspace(min(param[i]),max(param[i]),60),color='blue',label='string')
    plt.xlabel('Bins')
    plt.ylabel('Frequency')
    plt.title(titles[i])
