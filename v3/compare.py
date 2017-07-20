import numpy as np
import matplotlib.pyplot as plt
import sys

f1 = open(sys.argv[1], "r")
N = int(f1.readline())
T = int(f1.readline())

c1 = []

for j in range(1, N+1):
	cx1 = []	
	for i in range(1, N+1):	
		cx1.append(float(f1.readline()))
	c1.append(cx1)

f2 = open(sys.argv[2], "r")
N = int(f2.readline())
T = int(f2.readline())

c2 = []

for j in range(1, N+1):
	cx2 = []	
	for i in range(1, N+1):	
		cx2.append(float(f2.readline()))
	c2.append(cx2)


#h = np.array(c)
#print(np.linalg.norm(h) / np.sqrt(N))

h1 = np.array(c1)
h2 = np.array(c2)

plt.matshow(h1-h2)

plt.colorbar()

plt.show()

