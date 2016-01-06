import numpy as np
import itertools
import math
from matplotlib import pyplot as plt

def logpdf(data, step):
    xvalues = np.array(range(math.ceil(math.log10(max(data)) / math.log10(step)) + 1), dtype = float)
    xvalues *= math.log10(step)
    yvalues = xvalues * 0
    for value in data:
        if value >= 1:
            yvalues[math.ceil(math.log10(value) / math.log10(step))] += 1.0
        else:
            yvalues[0] += 1.0
    ysum = sum(yvalues)
    yvalues /= ysum
    print(xvalues, yvalues)
    return xvalues, yvalues


data = np.loadtxt("avalan")
flips = data[:,1] / 2
signed_flips = data[:,0]


print(max(flips))
print(math.log10(max(flips)))
x, y = logpdf(flips, 1.1)
sy = np.array([2,3,4,6,8,10], dtype = float)
sx = np.log10(sy)
sy /= 100
plt.plot(x, np.log10(y), '-ro', lw = 1)
#plt.plot(sx, sy,'bo')
plt.show()
