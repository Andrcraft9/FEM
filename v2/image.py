import numpy as np
import matplotlib.pyplot as plt
import sys

f = open(sys.argv[1], "r")

N = int(f.readline())
T = int(f.readline())

c = []

for j in range(1, N+1):
	cx = []	
	for i in range(1, N+1):	
		cx.append(float(f.readline()))
	c.append(cx)

plt.matshow(c)

plt.colorbar()

plt.show()

