import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

f = open(sys.argv[1], "r")
step = int(sys.argv[2])
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

im = plt.imshow(cT[step])

plt.colorbar()
plt.show()

