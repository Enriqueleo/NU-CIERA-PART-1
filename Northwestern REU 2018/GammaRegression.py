import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import random as r
import time
import scipy.stats as stats
s=time.time()
#create the x array from 1 to "" and the uncertainty of y
x=np.array(range(1,7051))/705
err=np.ones(len(x))/150



#This is the Gamma function
def gama(x,A,B,a):
    return A*(x**a)*np.exp(-B*x)
#This would be the "y" axis data values but this is randomized a little bit, 
#since gra is an array of random numbers. Ranx shifts the x values randomly
var=np.random.normal(0,err,len(x))
g=gama(x,100,5,6)+var
plt.figure()
plt.errorbar(x,g,yerr=err,fmt='o',color='blue')

cf,covr=spo.curve_fit(gama,x,g,[175,5,2],sigma=err)
gr=gama(x,cf[0],cf[1],cf[2])
chisq=sum(((g-gr)/err)**2)/(len(g)-3)
plt.plot(x,gr,'k')
erp=np.sqrt(np.diag(covr))
e=time.time()
print(co)
print(f'It took {e-s:1.2} seconds. ')
print(f'𝛘2={chisq:3.3}.')
print(f'A={cf[0]:4.4}±{erp[0]:1.3}, '\
f'β={cf[1]:4.4}±{erp[1]:1.2}, '\
f'and α={cf[2]:4.4}±{erp[2]:2.3}.')
