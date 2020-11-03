from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import math
import mpmath as mpm
import numpy as np
from scipy.linalg import expm

mpm.mp.dps = 30

n = 100
alpha = 32
k = 40
results_stronger = []
results_weaker = []
for delta in range(0, int(n/2)):
    alpha_over_k = mpm.mpf(alpha)/mpm.mpf(k)
    delta_over_n = mpm.mpf(delta)/mpm.mpf(n)
    # Get stronger result first
    p = mpm.mpf((int(n/2) - delta))/n
    t = mpm.mpf(alpha/k) - mpm.mpf(0.5) + mpm.mpf(delta/n)
    print(t)
    results_stronger.append(mpm.power(mpm.power(p/(p+t), p+t)*mpm.power((1-p)/(1-p-t), 1-p-t), k))
    # Get weaker results
    results_weaker.append(mpm.power(mpm.e, -2*mpm.power(alpha_over_k - mpm.mpf(0.5) + delta_over_n, 2)*k))

print(results_stronger)
plt.semilogy(range(0, int(n/2)), results_stronger, label="strong")
plt.semilogy(range(0, int(n/2)), results_weaker, label="weak")
plt.legend()
plt.show()