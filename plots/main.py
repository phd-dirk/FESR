import matplotlib.pyplot as plt
import numpy as np

# plt.plot([1,23,2,4])
# plt.ylabel('some numbers')

data = np.genfromtxt('../output/fits.dat', delimiter='\t', skip_header=1)
alpha = data[:,0]
aGGInv = data[:,2]
rhoVpA = data[:,4]
c8VpA = data[:, 6]
print(alpha)

plt.plot(alpha)
plt.plot(aGGInv)
plt.plot(rhoVpA)
plt.plot(c8VpA)

plt.show()
