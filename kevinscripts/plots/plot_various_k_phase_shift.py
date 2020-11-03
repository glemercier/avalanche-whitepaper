import matplotlib.pyplot as plt
import numpy as np
from sympy import Symbol, Rational, binomial, expand_func
import mpmath as mp
import math
import matplotlib as mpl
from multiprocessing import Pool
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

def test3Helper(t):
    n, k, alpha, c, b = t
    print(n, k, alpha, c, b)
    for i in range(0, 205):
        # print(i, ((0.5*c + i)*hyper(n, 0.5*c + i - 1, k, alpha)), ((0.5*c-i)*hyper(n, 0.5*c - i + b - 1, k, alpha)))
        if (((0.5*c-i)/c)*hyper(n, 0.5*c + i, k, alpha)) >= (((0.5*c+i)/c)*hyper(n, 0.5*c - i + b, k, alpha)):
            print(i)
            return i

def test3( _range, lab):
    r = []    
    with Pool(4) as p:
        r = p.map(test3Helper, _range)
    plt.plot([i[1] for i in _range], r, label=r"{}".format(lab))

def main():
    print(mpl.style.available)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    step = 2
    _range = [(1000, i, i, 800, 200) for i in range(1, 100+step, step)]
    test3(_range, "\alpha = 1")
    _range = [(1000, i, math.ceil(0.8*i), 800, 200) for i in range(1, 100+step, step)]
    test3(_range, "\alpha = 0.8")
    _range = [(1000, i, math.ceil(0.6*i), 800, 200) for i in range(1, 100+step, step)]
    test3(_range, "\alpha = 0.6")

    plt.xlabel(r'k', fontsize=60)
    plt.ylabel(r's_{ps}',fontsize=60)
    plt.legend()

    plt.show()


if __name__ == '__main__':
    main()
