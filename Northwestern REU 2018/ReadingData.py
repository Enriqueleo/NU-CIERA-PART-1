import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import glob
path1 = 'C:/Users/enriq/Desktop/Northwestern REU 2018/DATASET/sub DATA 1/*.txt'

#I read in the file of the Type Ia Supernova
#Fitting model equation
def f(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))
files1 = glob.glob(path1)

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
    bopt=-bopt*(1-0.14)/(1-0.02)
    bnd=np.array([(max(mag)-1.5,bopt-0.4,-1,-5,-5,4,t0-20),(max(mag)+1.5,bopt+0.4,1,5,5,70,t0+20)])
    pa,covar=spo.curve_fit(f,t,mag,[max(mag),bopt,-0.15,0.02,0.27,20,t0],sigma=err,bounds=bnd,method='trf')
    
    chisq.append(sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7))
    A.append(pa[0])
    b.append(pa[1])
    c.append(pa[2])
    d.append(pa[3])
    k.append(pa[4])
    sig.append(abs(pa[5]))
    tpeak.append(pa[6])
    x=np.linspace(min(t)-0.2,max(t)+1,100)
    plt.figure(figsize=(6,6))
    plt.plot(x,f(x,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),color='purple')
    plt.errorbar(t,mag,yerr=err,fmt='o',color='red')
    plt.ylabel('Magnitude',Fontsize='16',color='black')
    plt.xlabel('Epoch (days)',Fontsize='16',color='Black')
    plt.title("Light Curve",Fontsize='18',color='Black')
    plt.title(name)
    
A=np.array(A)
b=np.array(b)
c=np.array(c)
d=np.array(d)
k=np.array(k)
sig=np.array(sig)
tpeak=np.array(tpeak)
mean=[np.mean(A),np.mean(b),np.mean(c),np.median(d),np.median(k),np.median(sig),np.mean(tpeak)]
#This creates a histogram of the parameters from every data set that was 
#read in.
titles=['Amplitudes','b','c','d','k','σ','t0']
param=[A,b,c,d,k,sig,tpeak]
for i in range(len(titles)):
    fig=plt.figure(figsize=(9,8))
    plt.hist(param[i],bins=np.linspace(min(param[i]),max(param[i]),100),color='blue',label='string')
    plt.xlabel('Bins')
    plt.ylabel('Frequency')
    plt.title(titles[i])
