import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

f = open(sys.argv[1], "r")
N = int(f.readline())
T = int(f.readline())

cT = []
for k in range(1, T+1):
	c = []
	for j in range(1, N+1):
		cx = []	
		for i in range(1, N+1):	
			cx.append(float(f.readline()))
		c.append(cx)
	cT.append(c)


fig = plt.figure()
dt = 0
im = plt.imshow(cT[dt], animated=True)
#im_2 = plt.matshow(cT[dt], animated=True)


def updatefig(*args):
    global dt
    dt = (dt + 1) % T
    im.set_array(cT[dt])
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=True)

plt.colorbar()
plt.show()

