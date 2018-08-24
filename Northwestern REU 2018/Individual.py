import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import glob
path = 'C:/Users/enriq/Desktop/Northwestern REU 2018/DATASET/Interest not in use/*.txt'
files = glob.glob(path)
#I read in the file of the Type Ia Supernova
#Fitting model equation
def f(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))
for name in files:
    F=np.loadtxt(name)
    #The three col of information from the data
    t=F[:,0]
    mag=F[:,1]
    err=F[:,2]
    t0=t[np.where(mag==max(mag))[0]]
    bopt=(mag[len(t)-1]-mag[len(t)-3])/(t[len(t)-1]-t[len(t)-3])/max(mag)
    bopt=-bopt*(1-0.13)/(1-0.019)
    possig=t[len(t)-2]-t0
    lensig=possig/3
    possig=possig/4
    bnd=np.array([(max(mag)-1,bopt-0.2,-1,-0.5,0,0,t0-20),(max(mag)+1,bopt+0.2,0.5,1,4,possig+lensig,t0+20)])
    cf,covar=spo.curve_fit(f,t,mag,[max(mag),bopt,-0.13,0.019,0.3,possig,t0],sigma=err,bounds=bnd,method='trf')
    pa,covar=spo.curve_fit(f,t,mag,[max(mag),bopt,cf[2],0.02,0.3,possig,cf[6]],sigma=err,bounds=bnd,method='trf')
    plt.figure()
    plt.errorbar(t,mag,yerr=err,fmt='o')
    plt.title(name)
    x=np.linspace(min(t)-0.35,max(t)+1,100)
    plt.plot(x,f(x,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),color='purple')


    