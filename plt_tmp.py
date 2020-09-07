import numpy as np
import matplotlib.pyplot as plt
 

x = np.genfromtxt(r'tmp.dat')
plt.plot(x[:,0],x[:,2]-x[0,2],"ro--")
plt.show()
