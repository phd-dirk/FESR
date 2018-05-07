import matplotlib.pyplot as plt
import numpy as np

# plt.plot([1,23,2,4])
# plt.ylabel('some numbers')

data = np.genfromtxt('../output/fits.dat', delimiter='\t', skip_header=1)
alphas = data[:,0]
print(alphas)

plt.plot(alphas)

plt.show()
