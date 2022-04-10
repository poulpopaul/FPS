import matplotlib.pyplot as plt
import numpy as np

E = np.loadtxt("e.txt")
Ep = np.loadtxt("ep.txt")
Ec = np.loadtxt("ec.txt")
T = np.arange(0,len(E))
plt.plot(T,E,label = 'e')
plt.plot(T,Ec,label = 'ec')
plt.plot(T,Ep,label = 'ep')
#plt.ylim(-10,1E3)
#plt.yscale('log')
plt.legend()
plt.show()