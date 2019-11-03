import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text', usetex=True)

data1 = np.loadtxt("out/data0.33.csv", delimiter=",", skiprows=1)
data2 = np.loadtxt("out/data1.00.csv", delimiter=",", skiprows=1)
x1 = data1[:,0]
y1 = data1[:,1]
x2 = data2[:,0]
y2 = data2[:,1]

fig = plt.figure()
axes = plt.gca()
axes.set_ylim([0,120])
plt.plot(x1,y1,'-b',label=r'$\beta=0.33$')
plt.plot(x2,y2,'-r',label=r'$\beta=1.00$')
plt.legend()
plt.grid()
plt.xlabel("$x$ (m)")
plt.ylabel("$y$ (m)")
plt.title("Trajectories under different resistance laws")
fig.tight_layout()
fig.savefig('out/xy.png')
