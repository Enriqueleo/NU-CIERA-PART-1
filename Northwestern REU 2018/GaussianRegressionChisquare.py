import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import random as r
import time
s=time.time()
#create the x array from 1 to "" and the uncertainty of y
x=np.array(range(0,481))/24
err=np.ones(len(x))/5

#This is the Gaussian density function
def gauss(x,u,s,a):
    return (a/(s*np.sqrt(2*np.pi)))*np.exp((-1/2)*((x-u)/s)**2)
#This would be the "y" axis data values but this is randomized a little bit, 
#since gra is an array of random numbers
var=np.random.normal(0,1/5,len(x))
g=gauss(x,10,2,200)+var
#Plotting the data with it's uncertainties.
plt.figure()
plt.hist(var,bins=100)
plt.figure()
plt.errorbar(x,g,yerr=err,fmt='o',color='blue')
cf,covr=spo.curve_fit(gauss,x,g,[10,3,100],sigma=err)
gr=gauss(x,cf[0],cf[1],cf[2])
chisq=sum(((g-gr)/err)**2)/(len(g)-3)
erp=np.sqrt(np.diag(covr))
e=time.time()

print(covr)
print(f'It took {e-s:5.4} seconds.')
print(f'𝛘2={chisq:3.3}.')
print(f'μ={cf[0]:1.4}±{erp[0]:6.2}, '\
f'σ={cf[1]:1.3}±{erp[1]:6.2}, '\
f'and A={cf[2]:1.4}±{erp[2]:5.2}')
plt.plot(x,gr,'r')