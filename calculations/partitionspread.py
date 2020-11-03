import random
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
import math
from multiprocessing import Pool
from scipy.special import comb
from decimal import Decimal
import decimal
decimal.getcontext = 100000

# import functools
# import operator as op
# def ncr(n, r):
#     print(n, r)
#     r = min(r, n-r)
#     if r == 0: return 1
#     numer = functools.reduce(op.mul, range(n, n-r, -1))
#     denom = functools.reduce(op.mul, range(1, r+1))
#     return numer//denom

def ncr(n, r):
    return comb(n, r, exact=True)

def hyper(N, succ, sample, x):
    succi = int(succ)
    outcome = Decimal(0)
    for i in range(x, sample+1):
        # print(":::", succi, i, N-succi, sample-i, sample, N)
        top = ncr(succi, i) * ncr(N - succi, sample - i)
        bottom = ncr(N, sample)
        outcome += (Decimal(top)/Decimal(bottom))
    return outcome
    # rv = hypergeom(N, succi, sample)
    # return decimal.Decimal(1.0) - decimal.Decimal(rv.cdf(x)) + decimal.Decimal(rv.pmf(x))

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
    T1sum = Decimal(0)
    T2sum = Decimal(0)
    # Time 0, base case
    gamma1 = Decimal(Z1/C)*hyper(N, C1 + B, k, alpha)
    gamma2 = Decimal(Z2/C)*hyper(N, C1, k, alpha)
    T1sum += (Decimal(gamma1) + Decimal(gamma2))
    xi1 = Decimal(Z1/C)*hyper(N, C2, k, alpha)
    xi2 = Decimal(Z2/C)*hyper(N, C2 + B, k, alpha)
    T2sum += (Decimal(xi1) + Decimal(xi2))
    gammas = [Decimal(gamma2)]
    xis = [Decimal(xi2)]
    # other cases
    for i in range(1, C):
        # print(i, T1sum, T2sum)
        try:
            assert Decimal(C) > T1sum + T2sum
        except Exception as e:
            print(T1sum, T2sum)
            raise e
        gamma1 = Decimal(Z1/C)*hyper(N, C1 + B + sum(gammas) - sum(xis), k, alpha)
        # print(Decimal(Z1/C), hyper(N, C1 + B + sum(gammas) - sum(xis), k, alpha))
        gamma2 = Decimal(Z2/C)*hyper(N, C1 + sum(gammas) - sum(xis), k, alpha)
        # print(Decimal(Z2/C), hyper(N, C1 + sum(gammas) - sum(xis), k, alpha))
        xi1 = Decimal(Z1/C)*hyper(N, C2 + sum(xis) - sum(xis), k, alpha)
        xi2 = Decimal(Z2/C)*hyper(N, C2 + B + sum(xis) - sum(xis), k, alpha)
        # print(gamma1, gamma2, xi1, xi2)
        if math.isnan(gamma1):
            gamma1 = 0.0
        if math.isnan(gamma2):
            gamma2 = 0.0
        if math.isnan(xi1):
            xi1 = 0.0
        if math.isnan(xi2):
            xi2 = 0.0
        T1sum += (Decimal(gamma1) + Decimal(gamma2))
        T2sum += (Decimal(xi1) + Decimal(xi2))
        gammas.append(Decimal(gamma2))
        xis.append(Decimal(xi2))
    return (t[4], T1sum/C, T2sum/C)

def main():
    k = 40
    alpha = 24
    N = 1000
    B = 200
    with Pool(8) as p:
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(0, N - B + 10, 10)])
    # r = [exp(200)]
    x = [i[0] for i in r]
    y1 = [0 if math.isnan(i[1]) else i[1] for i in r]
    y2 = [0 if math.isnan(i[2]) else i[2] for i in r]
    plt.plot(x, y1, label="C1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))
    plt.plot(x, y2, label="C2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))

    k = 40
    alpha = 28
    N = 1000
    B = 200
    with Pool(8) as p:
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(0, N - B + 10, 10)])
    # r = [exp(200)]
    x = [i[0] for i in r]
    y1 = [0 if math.isnan(i[1]) else i[1] for i in r]
    y2 = [0 if math.isnan(i[2]) else i[2] for i in r]
    plt.plot(x, y1, label="C1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))
    plt.plot(x, y2, label="C2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))

    k = 40
    alpha = 32
    N = 1000
    B = 200
    with Pool(8) as p:
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(0, N - B + 10, 10)])
    # r = [exp(200)]
    x = [i[0] for i in r]
    y1 = [0 if math.isnan(i[1]) else i[1] for i in r]
    y2 = [0 if math.isnan(i[2]) else i[2] for i in r]
    plt.plot(x, y1, label="C1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))
    plt.plot(x, y2, label="C2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))

    k = 40
    alpha = 36
    N = 1000
    B = 200
    with Pool(8) as p:
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(0, N - B + 10, 10)])
    # r = [exp(200)]
    x = [i[0] for i in r]
    y1 = [0 if math.isnan(i[1]) else i[1] for i in r]
    y2 = [0 if math.isnan(i[2]) else i[2] for i in r]
    plt.plot(x, y1, label="C1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))
    plt.plot(x, y2, label="C2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B))

    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()