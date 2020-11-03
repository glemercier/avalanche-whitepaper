import numpy as np
from multiprocessing import Pool
import math
from sympy import binomial, zeros
import mpmath as mp
mp.mp.dps = 20
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.style.use('seaborn-white')
import numpy as np
import pickle

def hyper(N, suc, sample, x):
    total = mp.mp.mpf(0)
    bottom = binomial(N, sample)
    for i in range(x, sample+1):
        top = binomial(suc, i) * binomial(N-suc, sample-i)
        total += (top/bottom)
    return total

def birthDeathProbability(n, k, alpha, c, b):
    P = zeros(c/2, c/2)
    for i in range(1, int(0.5*c)-1):
        print(i)
        birth = ((0.5*c-i)/c)*hyper(n, 0.5*c + i, k, alpha)
        death = ((0.5*c+i)/c)*hyper(n, 0.5*c - i + b, k, alpha)
        P[i, i-1] = death
        P[i, i+1] = birth
        P[i, i] = mp.mp.mpf(1.0) - birth - death
    # Half blue half red
    # Choose blue, get red
    P[0, 1] = (0.5)*hyper(n, c/2, k, alpha)
    # Choose blue, get nothing
    P[0, 0] = mp.mp.mpf(1.0) - (0.5)*hyper(n, c/2, k, alpha)
    # No blue, all red
    # Choose red, get blue
    P[int(0.5*c)-1, int(0.5*c)-2] = hyper(n, b, k, alpha)
    # Choose red, get nothing
    P[int(0.5*c)-1, int(0.5*c)-1] = mp.mp.mpf(1.0) - hyper(n, b, k, alpha)
    powerraise = 50000
    P = P**powerraise
    with open("./P_{}_{}_{}_{}_{}_{}.out".format(n,k,alpha,c,b,powerraise), 'wb+') as f:
        pickle.dump(P, f)
    # plt.plot(range(0, int(c/2)), P.col(20))
    P = np.array(P).astype(np.float64)    

    x = range(0, int(c/2))
    y = range(0, int(c/2))

    X, Y = np.meshgrid(x, y)
    Z = P

    # plt.imshow(Z, extent=[0, int(c/2), 0, int(c/2)], origin='lower', cmap=plt.cm.cool, alpha=0.5, norm=LogNorm())
    plt.matshow(P, cmap=plt.cm.gnuplot2)
    # contours = plt.contour(X, Y, P, cmap=plt.cm.cool)
    # plt.clabel(contours, inline=True, fontsize=12)
    plt.colorbar()

    plt.show()

def main():
    if True:
        i = 10
        birthDeathProbability(10*i, 20, 16, 8*i, 2*i)
        exit()

if __name__ == '__main__':
    main()