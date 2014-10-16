#!/usr/bin/python
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

np.loadtxt("vanderpol.txt")
data = np.loadtxt("vanderpol.txt")
xsol = data[:,0]
ysol = data[:,1]
dt   = 0.2
nt   = np.size(xsol)
t    = np.linspace(0, nt*dt, nt)

plt.figure(figsize=(8,8))
plt.plot(t, xsol, color='red', label='x')
plt.plot(t, ysol, color='blue', label='y')
plt.xlim(0, t[nt-1])
plt.legend()
plt.show()