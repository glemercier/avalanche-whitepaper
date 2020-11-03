import random
import math
import itertools
from scipy.stats import hypergeom
from functools import reduce
import operator

def hyper(N, S, k, a):
    rv = hypergeom(N, S, k)
    return 1.0 - rv.cdf(a) + rv.pmf(a)

def prod(iterable):
    if not len(iterable):
        return 1
    return reduce(operator.mul, iterable, 1)

def t1(births, deaths):
    sigma = 0
    for k in range(1, len(deaths)):
        top = prod(births[0:k])
        bottom = prod(deaths[0:k+1])
        sigma += (top/bottom)
    return (1/deaths[0])


# def extinction_time(births, deaths):
#     stays = [1-(births[i]+deaths[i]) for i in range(len(births))]
#     print(births, deaths, stays)
#     tm = t1(births, deaths)

def An(n, p, q):
    sigma = 0
    for i in range(1, n-1):
        print(p[1:i+1])
        print(q[1:i+1])
        print("-----------")
        top = prod(q[1:i+1])
        bottom = prod(p[1:i+1])
        sigma += top/bottom
    return sigma

def mn1(n, p, q):
    a = 1/An(n, p, q)
    top_sigma = 0
    for y in range(1, n-1):
        bottom_sigma = 0
        for z in range(1, y+1):
            top = prod(q[z+1:y+1])
            bottom = prod(p[z:y+1])
            print("y {}, z {}, top {}, bottom {}".format(y, z, q[z+1:y+1], p[z:y+1]))
            bottom_sigma += (top/bottom)
        top_sigma += bottom_sigma
    print(a * top_sigma)
    return a * top_sigma

def mn(x, p, q):
    _mn1 = mn1(len(p), p, q)
    

# Birth death discrete chain with absorbing endpoints
def vn(x, p, q):
    f = 1/An(len(p), p, q)
    sigma = 0
    for i in range(x, len(p)-2):
        top = prod(q[1:i+1])
        bottom = prod(p[1:i+1])
        sigma += (top/bottom)
    return(f*sigma)

print(mn(2, [0, 1/2, 1/2, 1/2, 0], [0, 1/2, 1/2, 1/2, 0]))