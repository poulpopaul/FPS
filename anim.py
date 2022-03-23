import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

X = np.loadtxt(r"x.txt")
Y = np.loadtxt(r"y.txt")

n = X[0]


ims = []
fig = plt.figure()
plt.xlim(0, 32)
plt.ylim(0, 32)
for k in range(100):
    a = 1 +k*int(n)
    b = int(n)+ k*int(n)
    ims.append(plt.plot(X[a:b],Y[a:b],linestyle = ' ',marker = 'o',color = 'r',alpha = 0.5,markersize = 2))

plt.show()

ani = animation.ArtistAnimation(fig, ims,interval = 1,repeat = True)
ani.save("ani.gif",fps = 1000,writer = 'Pillow')
#plt.show(block = False)