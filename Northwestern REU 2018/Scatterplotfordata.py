import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import pandas as pd
import glob
F=np.loadtxt\
('DATASET/Look at later/SN2004gcopt_gFilter_.txt')
def f(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))
t=F[:,0]
mag=F[:,1]
err=F[:,2]
plt.scatter(t,mag)
t0=t[np.where(mag==max(mag))[0]]
bopt=(mag[len(t)-1]-mag[len(t)-3])/(t[len(t)-1]-t[len(t)-3])/max(mag)
pa,covar=spo.curve_fit(f,t,mag,[max(mag),0.02,-0.0005,-0.15,0.4,30,t0],sigma=err)

plt.plot(t,f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]),color='purple')
chisq=sum(((mag-f(t,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6]))/err)**2)/(len(t)-7)
print(chisq)
plt.savefig()
