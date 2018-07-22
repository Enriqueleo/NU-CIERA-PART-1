import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import glob
import errno
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

fig=plt.figure(figsize=(13,12))
plt.hist(A,bins=np.linspace(-20,-13,29),color='red',label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('Amplitudes')
fig.savefig('Amp.pdf')

fig=plt.figure(figsize=(13,12))
plt.hist(b,bins=np.linspace(min(b),max(b),40),color='black',label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('b')
fig.savefig('b.pdf')

fig=plt.figure(figsize=(13,12))
plt.hist(c,bins=np.linspace(min(c),max(c),50),color='Green',label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('c')
fig.savefig('c.pdf')

fig=plt.figure(figsize=(13,12))
plt.hist(d,bins=np.linspace(min(d),max(d),50),color='red',label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('d')
fig.savefig('d.pdf')

fig=plt.figure(figsize=(13,12))
plt.hist(k,bins=np.linspace(min(k),max(k),50),color='Navy',stacked=True,label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('k')
fig.savefig('k.pdf')

fig=plt.figure(figsize=(13,12))
plt.hist(sig,bins=np.linspace(min(sig),max(sig),50),color='Purple',label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('σ')
fig.savefig('sigma.pdf')

fig=plt.figure(figsize=(13,12))
plt.hist(tpeak,bins=np.linspace(min(tpeak),max(tpeak),50),color='Pink',label='string')
plt.xlabel('Bins')
plt.ylabel('Frequency')
plt.title('t0')
fig.savefig('tatmax.pdf')