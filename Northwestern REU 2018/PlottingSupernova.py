#Set up plotting libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import glob
#Set up astropy 
from astropy.coordinates import SkyCoord, Distance
from astropy import units as u
import re

#I read in the file of the Type Ia Supernova
#Fitting model equation
def LC(xraw,A,b,c,d,k,s,t):
    x=xraw-t
    return ((A-A*d)/(1+c))*((1-b*x+c*np.exp(-(x/s)**2))/(1-d*np.exp(-k*x)))

deltaM=[]
Mmax=[]

RA=[]
DEC=[]
#So distance is porportional to resshift. The larger the redshift
#the further the distance. Inverting the redshift means small redshifts( closer)
#bigger the data point.
Zinverse=[]
X=[]
Y=[]
info=[]

path='More DATA/Nor/*.dat'
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
    cf,covar=spo.curve_fit(LC,t,mag,[min(mag),bopt,-0.13,0.018,0.3,possig,t0],sigma=err,bounds=bnd,method='trf')
    bnd=np.array([(min(mag)-1,cf[1]-0.1,cf[2]-0.1,cf[3]-0.5,0,0,t0-7),(min(mag)+1,cf[1]+0.2,cf[2]+0.1,1,5,possig+lensig,t0+7)])
    cf,covar=spo.curve_fit(LC,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd,method='trf')
    possig=cf[5]
    bnd=np.array([(min(mag)-1,cf[1]-0.05,cf[2]-0.1,cf[3]-0.25,0,0,t0-2),(min(mag)+1,cf[1]+0.3,cf[2]+0.1,cf[3]+0.05,1.5,possig+12,t0+2)])
    pa,covar=spo.curve_fit(LC,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd)
    
    DM=LC(pa[6]+15,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])-LC(pa[6],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])
    deltaM.append(DM)
    Mmax.append(-21.726+2.698*DM)
    
    item=[]
    for line in open(name, 'r'):
        item.append(line.rstrip())
    info.append(item[2])
    num = np.loadtxt(re.findall(r"[-+]?\d*\.\d+|\d+",item[2]))
    RA.append(15*(num[1]+num[2]/60+num[3]/3600-12))
    X.append(np.pi*(num[1]+num[2]/60+num[3]/3600)/12)
    DEC.append(num[4]+num[5]/60+num[5]/3600)
    Y.append(-np.cos(np.pi*(90+(num[4]+num[5]/60+num[5]/3600))/180))
    Zinverse.append(6/num[0])

path='More DATA/Sur/*.dat'
file2=glob.glob(path)

for name in file2:
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
    cf,covar=spo.curve_fit(LC,t,mag,[min(mag),bopt,-0.13,0.018,0.3,possig,t0],sigma=err,bounds=bnd,method='trf')
    bnd=np.array([(min(mag)-1,cf[1]-0.1,cf[2]-0.1,cf[3]-0.5,0,0,t0-7),(min(mag)+1,cf[1]+0.2,cf[2]+0.1,1,5,possig+lensig,t0+7)])
    cf,covar=spo.curve_fit(LC,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd,method='trf')
    possig=cf[5]
    bnd=np.array([(min(mag)-1,cf[1]-0.05,cf[2]-0.1,cf[3]-0.25,0,0,t0-2),(min(mag)+1,cf[1]+0.3,cf[2]+0.1,cf[3]+0.05,1.5,possig+12,t0+2)])
    pa,covar=spo.curve_fit(LC,t,mag,[min(mag),cf[1],cf[2],cf[3],cf[4],possig,t0],sigma=err,bounds=bnd)
    
    DM=LC(pa[6]+15,pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])-LC(pa[6],pa[0],pa[1],pa[2],pa[3],pa[4],pa[5],pa[6])
    deltaM.append(DM)
    Mmax.append(-21.726+2.698*DM)
    
    item=[]
    for line in open(name, 'r'):
        item.append(line.rstrip())
    info.append(item[2])
    num = np.loadtxt(re.findall(r"[-+]?\d*\.\d+|\d+",item[2]))
    RA.append(15*(num[1]+num[2]/60+num[3]/3600-12))
    X.append(np.pi*(num[1]+num[2]/60+num[3]/3600)/12)
    DEC.append(-(num[4]+num[5]/60+num[5]/3600))
    Y.append(-np.cos(np.pi*(90-(num[4]+num[5]/60+num[5]/3600))/180))
    Zinverse.append(6/num[0])

f, ax = plt.subplots(subplot_kw={'projection': "mollweide"},figsize=(10,8))
ax.grid(True)
Coords=SkyCoord(RA,DEC,unit=(u.degree, u.degree))
ax.scatter(Coords.ra.wrap_at(180.*u.degree).radian,Coords.dec.radian,c=Mmax,cmap='cool',s=Zinverse,linewidths=1.25,edgecolors='black')
ax.set_xlabel("RA",fontsize=20)
ax.set_ylabel("DEC",fontsize=20)
plt.title('Type Ia Supernova Locations',fontsize=25)
plt.savefig("Supernova location.png")
plt.figure(figsize=(10,8))
plt.scatter(X,Y,cmap='cool',c=Mmax,s=Zinverse,linewidths=1.25,edgecolors='black')
plt.xlim([-.1,2.05*np.pi])
plt.ylim([-1,1])
plt.ylabel('-Cos(DEC+90)')
plt.xlabel('RA (radians)')
plt.title('Distribution of supernova',fontsize=15)
plt.figure()
plt.hist(X,bins=np.linspace(0,2*np.pi,20))
plt.title('RA')
plt.xlabel('Radians')
plt.figure()
plt.hist(Y,bins=np.linspace(-1,1,21))
plt.title("-Cos(DEC+90)")
