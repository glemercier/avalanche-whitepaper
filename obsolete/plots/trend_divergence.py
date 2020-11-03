import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from sympy import binomial, zeros, pprint
import mpmath as mp
import math
mp.mp.dps = 20

def hyper(N, suc, sample, x):
    total = mp.mp.mpf(0)
    bottom = binomial(N, sample)
    for i in range(x, sample+1):
        top = binomial(suc, i) * binomial(N-suc, sample-i)
        total += (top/bottom)
    return total

def main():
    n = 1000
    b = 200
    k = 10
    alpha = 8
    _range = range(0, 1001)
    plot1 = [hyper(1000, i, k, alpha) for i in _range]
    plot2 = [hyper(1000, 1000-i, k, alpha) for i in _range]
    plt.semilogy(_range, plot1, color="red", label="Probability of \mathtt{R} successful sample")
    plt.semilogy(_range, plot2, color="blue", label="Probability of \mathtt{B} successful sample")
    plt.xlabel("Number of \mathtt{R}")
    plt.ylabel("Probability of a success")
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()