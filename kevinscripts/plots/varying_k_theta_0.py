import matplotlib.pyplot as plt
import numpy as np
from sympy import Symbol, Rational, binomial, expand_func
import mpmath as mp
import math
import matplotlib as mpl
# mpl.style.use("seaborn-dark-palette")

params = {'legend.fontsize': 40,
          'figure.figsize': (15, 5),
         'axes.labelsize': 60,
         'axes.titlesize': 60,
         'xtick.labelsize': 60,
         'ytick.labelsize': 60}
plt.rcParams.update(params)

mp.mp.dps = 20

def hyper(N, suc, sample, x):
    total = mp.mp.mpf(0)
    bottom = binomial(N, sample)
    for i in range(x, sample+1):
        top = binomial(suc, i) * binomial(N-suc, sample-i)
        total += (top/bottom)
    return total

print(mpl.style.available)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
for k in [1, 5, 10, 20, 40, 80, 160]:
    print(k)
    alpha = math.ceil(0.8*k)
    plt.semilogy(range(400, 800), [hyper(1000, i, k, alpha) for i in range(400, 800)], label=r"k = {}".format(k))

plt.xlabel(r's_i', fontsize=60)
plt.ylabel(r'p(s_i)',fontsize=60)
plt.legend()

plt.show()
