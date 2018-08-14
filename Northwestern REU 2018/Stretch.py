import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import glob
import pandas as pa
import re

Z=[]
SSF=[]
info=[]

deltaM=[]
Mmax=[]
A=[]
b=[]
c=[]
d=[]
k=[]
sig=[]
tpeak=[]
ESOP=[]
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
    lensig=possig/3.333
    possig=possig/4
    bnd=np.array([(min(mag)-1,bopt-0.15,-0.5,-0.3,0,0,t0-15),(min(mag)+1,bopt+0.15,0.5,0.5,2,possig+lensig,t0+15)])
    cf,covar=spo.curve_fit(f,t,mag,[min(mag),bopt,-0.13,0.018,0.3,possig,t0],sigma=err,bounds=bnd,method='trf')
    bnd=np.array([(min(mag)-1,cf[1]-0.1,cf[2]-0.1,cf[3]-0.5,0,0,t0-7),(min(mag)+1,cf[1]+0.2,cf[2]+0.1,1,5,possig+lensig,t0+7)])
    cf,covar=spo.curve_fit(f,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd,method='trf')
    possig=cf[5]
    bnd=np.array([(min(mag)-1,cf[1]-0.05,cf[2]-0.1,cf[3]-0.25,0,0,t0-2),(min(mag)+1,cf[1]+0.3,cf[2]+0.1,cf[3]+0.05,1.5,possig+12,t0+2)])
    pa,covar=spo.curve_fit(f,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd)
    A.append(pa[0])
    b.append(pa[1])
    c.append(pa[2])
    d.append(pa[3])
    k.append(pa[4])
    sig.append(pa[5])
    tpeak.append(pa[6])
    ESOP.append(pa)
    item=[]
    for line in open(name, 'r'):
        item.append(line.rstrip())
    info.append(item[2])
    num = np.loadtxt(re.findall(r"[-+]?\d*\.\d+|\d+",item[2]))
    Z.append(num[0])
    
    deltaM.append(f(pa[6]+15,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])-f(pa[6],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))
    Mmax.append(-21.726+2.698*(f(pa[6]+15,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])-f(pa[6],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])))

MeanA=np.mean(np.array(A))
Meanb=np.mean(np.array(b))
Meanc=np.mean(np.array(c))
Meand=np.mean(np.array(d))
Meank=np.mean(np.array(k))
Meansig=np.mean(np.array(sig))


m15bar=f(15,MeanA,Meanb,Meanc,Meand,Meank,Meansig,0)
x=np.linspace(0,100,3502)
for i in range(len(ESOP)):
    QQ=np.where(m15bar<=f(x/(1+Z[i]),ESOP[i][0],ESOP[i][1],ESOP[i][2],ESOP[i][3],ESOP[i][4],ESOP[i][5],0)+MeanA-ESOP[i][0])[0][0]
    SSF.append(x[QQ]/(15+15*Z[i]))

j=0
legen=[]
plt.figure(figsize=(16,14))
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
    leg=name.replace('More DATA/Graph\\','')
    legen.append(leg.replace('+nir_photo.dat',''))
    t=np.array(t)
    mag=np.array(mag)
    err=np.array(err)
    te=(t-ESOP[j][6])/((1+Z[j])*SSF[j])
    plt.scatter(te,mag+MeanA-ESOP[j][0],s=50)
    plt.xlabel('Days',fontsize='15')
    plt.ylabel('B + offset',fontsize='15')
    j+=1
plt.gca().invert_yaxis()
plt.xlim([-20,100])
plt.ylim([20,15])
plt.title('Standardized',fontsize='25')
plt.legend(legen)
