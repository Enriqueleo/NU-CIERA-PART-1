import matplotlib.pyplot as plt
import glob
import pandas as pa
import numpy as np
path='More DATA/Graph/*.dat'
file=glob.glob(path)
for name in file:
    t=[]
    mag=[]
    err=[]
    F=np.loadtxt(name)
    myarray=pa.read_csv(name)
    for i in range(len(F)):
        if F[i,4]<=99:
            t.append(F[i,0])
            mag.append(-F[i,3])
            err.append(F[i,4])
    plt.figure()
    plt.scatter(t,mag)
    plt.title(name)
item=[]
for line in open(name, 'r'):
    item.append(line.rstrip())
import re

