import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
print("The GIF animation is being created...\n")
X = np.loadtxt(r"x.txt")
Y = np.loadtxt(r"y.txt")
T = np.loadtxt(r"t.txt")
C = [1]
for t in T:
    if t==T[0]:
        C.append('r')
    else:
        C.append('b')

n = X[0]
frame = int(len(X)/n)
x_max = max(X[1:int(n)])
y_max = max(Y[1:int(n)])


fig, ax = plt.subplots(figsize=(6,6))
def update(i):
    ax.clear()
    ax.set_xlim(0,x_max)
    ax.set_ylim(0,y_max)
    ax.set_xticks([])
    ax.set_yticks([])
    a = 1 +i*int(n)
    b = int(n)+ i*int(n) +1
    ax.scatter(X[a:b],Y[a:b],marker = 'o',c = C[a:b],alpha = 0.6,s = 1)


ani = animation.FuncAnimation(fig, update, frames=frame,interval = frame/4)
writergif = animation.PillowWriter(fps = 30) 
ani.save("ani.gif",writer=writergif)
print("ani.gif has been saved !\n")