import random
import math
from multiprocessing import Pool

import numpy as npy
import mpmath as mpmath
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from scipy.special import comb

mpmath.mp.dps = 200

def ncr(n, r):
    return comb(n, r, exact=True)

def hyper(N, succ, sample, x):
    outcome = mpmath.mpf(0)
    for i in range(x, sample+1):
        # print(":::", succi, i, N-succi, sample-i, sample, N)
        top = ncr(succ, i) * ncr(N - succ, sample - i)
        bottom = ncr(N, sample)
        outcome += (mpmath.mpf(top)/mpmath.mpf(bottom))
    return outcome

def exp(t):
    k = t[0]
    alpha = t[1]
    N = t[2]
    B = t[3]
    C = N - B
    C1 = t[4]
    C2 = C - C1
    Z1 = C1
    Z2 = C2
    T1sum = mpmath.mpf(0)
    T2sum = mpmath.mpf(0)
    # Generate matrix
    npy.matrix
    return (t[4], T1sum, T2sum)

def main():
    k = 40
    alpha = 32
    N = 1000
    B = 200
    step = 100
    with Pool(8) as p:
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(0, int((N - B)/2) + step, step)])
    # r = [exp(200)]
    x = [i[0] for i in r]
    y1 = [0 if math.isnan(i[1]) else i[1] for i in r]
    y2 = [0 if math.isnan(i[2]) else i[2] for i in r]
    plt.plot(x, y1, label="C1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))
    plt.plot(x, y2, label="C2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))
    print(min(y1), max(y1))
    print(min(y2), max(y2))
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()