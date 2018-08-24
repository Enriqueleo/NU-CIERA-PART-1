import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import glob
import pandas as pa

A=[]
b=[]
c=[]
d=[]
k=[]
sig=[]
tpeak=[]
chisq=[]
EachSetOfPara=[]

deltaM=[]
Mmax=[]

#I read in the file of the Type Ia Supernova
#Fitting model equation
def f(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))

path='More DATA/Graph/*.dat'
file=glob.glob(path)
for name in file:
    t=[]
    mag=[]
    err=[]
    F=np.loadtxt(name)    
    for i in range(len(F)):
        if F[i,4]<99:
            t.append(F[i,0])
            mag.append(F[i,3])
            err.append(F[i,4])
    t=np.array(t)
    mag=np.array(mag)
    err=np.array(err)
    t0=t[np.where(mag==min(mag))[0]]
    bopt=(mag[len(t)-2]-mag[len(t)-4])/(t[len(t)-2]-t[len(t)-4])/min(mag)
    bopt=-bopt*(1-0.13)/(1-0.018)
    possig=t[len(t)-3]-t0
    lensig=possig/3.3333
    possig=possig/4
    bnd=np.array([(min(mag)-1,bopt-0.15,-0.5,-0.3,0,0,t0-15),(min(mag)+1,bopt+0.15,0.5,0.5,2,possig+lensig,t0+15)])
    cf,covar=spo.curve_fit(f,t,mag,[min(mag),bopt,-0.13,0.018,0.3,possig,t0],sigma=err,bounds=bnd,method='trf')
    bnd=np.array([(min(mag)-1,cf[1]-0.1,cf[2]-0.1,cf[3]-0.5,0,0,t0-7),(min(mag)+1,cf[1]+0.2,cf[2]+0.1,1,5,possig+lensig,t0+7)])
    cf,covar=spo.curve_fit(f,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd,method='trf')
    possig=cf[5]
    bnd=np.array([(min(mag)-1,cf[1]-0.05,cf[2]-0.1,cf[3]-0.25,0,0,t0-2),(min(mag)+1,cf[1]+0.3,cf[2]+0.1,cf[3]+0.05,1.5,possig+12,t0+2)])
    pa,covar=spo.curve_fit(f,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd)
    
    x=np.linspace(min(t)-0.5,max(t)+5,1501)
    plt.figure(figsize=(12,8))
    plt.errorbar(t-pa[6],mag,yerr=err,fmt='o',capsize=15,color='red')
    plt.scatter(t-pa[6],mag,marker='o',s=50,c='r')
    plt.plot(x-pa[6],f(x,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),color='purple')
    plt.xlabel('Days',fontsize=18)
    plt.ylabel('B',fontsize=18)
    plt.gca().invert_yaxis()
    chisq.append(sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7))
    A.append(pa[0])
    b.append(pa[1])
    c.append(pa[2])
    d.append(pa[3])
    k.append(pa[4])
    sig.append(pa[5])
    tpeak.append(pa[6])
    EachSetOfPara.append(pa)
    deltaM.append(f(pa[6]+15,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])-f(pa[6],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))
    Mmax.append(-21.726+2.698*(f(pa[6]+15,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])-f(pa[6],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])))
    leg=name.replace('More DATA/Graph\\','')
    leg=leg.replace('+nir_photo.dat','')
    plt.title(leg,fontsize=20)
plt.savefig(leg)
titles=['b','c','d','k','σ','Δm15(B)']
param=[b,c,d,k,sig,deltaM]
for i in range(len(titles)):
    plt.figure(figsize=(10,5))
    plt.hist(param[i],bins=np.linspace(min(param[i]),max(param[i]),51),color='red',label='string')
    plt.xlabel('Bins')
    plt.ylabel('Frequency')
    plt.title(titles[i],fontsize=18)

plt.figure(figsize=(12,8))
plt.plot(deltaM,Mmax,color='pink')
plt.scatter(deltaM,Mmax,s=100,c='purple')
plt.ylabel('Mmax(B)',fontsize=15)
plt.xlabel('Δm15(B)',fontsize=15)
plt.title("Phillips' Relationship",fontsize=18)
plt.gca().invert_yaxis()
plt.savefig("Phillips Relationship.png")