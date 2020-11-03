import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import hypergeom

def hyper(N, suc, sample, x):                                                                                                   
    rv = hypergeom(N, suc, sample)                                                                                              
    return 1.0 - rv.cdf(x) + rv.pmf(x)

plt.semilogx([10**i for i in range(2, 10)], [hyper(10**i, int((10**i)*0.8), 5, 4)/hyper(100, 80, 5, 4) for i in range(2, 10)])
plt.xlabel("Network size")
plt.ylabel("Ratio")
plt.title("Ratios of the probabilities of successes over various network sizes, with fixed k = 10, alpha = 0.8, and x = 0.8c.")
plt.show()