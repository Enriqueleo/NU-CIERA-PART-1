import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import random as r

x=np.array(range(-10,71))/10
err=np.ones(len(x))
y=2*x+10+np.random.normal(0,err,len(x))
def linModel(u,m,c):
    return m*u+c


plt.figure()
plt.errorbar(x,y,yerr=err,fmt='o',color='blue')

coef,covr=spo.curve_fit(linModel,x,y,sigma=err)
yr=linModel(x,coef[0],coef[1])

errp=np.sqrt(np.diag(covr))
plt.plot(x,yr,'r')
chisq = sum(((yr-y)/err)**2)/(len(y)-2)
print(f'The chi square is {chisq :4.4}, '\
f'the slope is {coef[0]:4.4}, +/- {errp[0]:4.2} '\
f'and the y-intercept is {coef[1]:4.4}+/-{errp[1]:4.2}.')