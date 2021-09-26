import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro')
xdata2, ydata2 = [], []
ln2, = plt.plot([], [], 'b-')

N=1000
dim=2
f = open("coordinates.txt", "r")
coor_x=[]
coor_y=[]
coor_z=[]

for line in f:
    line=line.split('\t')
    coor_x.append(float(line[0]))
    coor_y.append(float(line[1]))
    if dim == 3:
        coor_z.append(float(line[2]))

def init():
    ax.set_xlim(-100, 100)
    ax.set_ylim(-100, 100)
    return ln,

def update(frame):
    xdata.append(coor_x[frame])
    ydata.append(coor_y[frame])
    ln.set_data(xdata, ydata)
    if frame>0:
        xdata2.extend(np.linspace(coor_x[frame-1],coor_x[frame],10))
        ydata2.extend(np.linspace(coor_y[frame-1],coor_y[frame],10))
    else:
        xdata2.append(coor_x[frame])
        ydata2.append(coor_y[frame])
    ln2.set_data(xdata2, ydata2)
    return ln, ln2,


ani = FuncAnimation(fig, update, frames=np.arange(N), init_func=init, blit=True)
plt.show()