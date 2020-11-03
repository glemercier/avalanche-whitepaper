import multiprocessing
from decimal import Decimal
import decimal
decimal.getcontext = 1000
import matplotlib.pyplot as plt
import math
import numpy as np

def nck(n, k):
    return Decimal(math.factorial(n))/(Decimal(math.factorial(k))*Decimal(math.factorial(n-k)))

def hyperpmf(population, successes, draws, observed_successes):
    ret = nck(successes, observed_successes)*(nck(population-successes, draws-observed_successes)/nck(population, draws))
    return ret

def hyperormore(tup):
    total = Decimal('0')
    N = tup[0]
    K = tup[1]
    n = tup[2]
    k = tup[3]
    for i in range(k, n+1):
        total += hyperpmf(N, K, n, i)
    return total

# with multiprocessing.Pool(20) as p:
#     ranges = range(5000, 8200, 200)
#     y = p.map(hyperormore, [(10000, i, 100, 80) for i in ranges])

# plt.plot(ranges, y)
# plt.show()

from scipy.stats import hypergeom

def hyperormore_fast(pop, succ, sample, threshold):
    rv = hypergeom(pop, succ, sample)
    return 1 - rv.cdf(threshold) + rv.pmf(threshold)

import pykov

'''
Attempting to instantiate the Markov process that describes an attack where 
the Byzantine nodes partition the network into some H1 and H2 sets, and then 
attempt to constantly try to prefer one trasanaction over the other. 
'''
import random
data = []
trials = 10000
ksample = 50
alpha = 40
total_nodes = 100000
total_byz = 20000
ranges = range(40000, 80000, 50)
for i in ranges:
    total_h1 = i
    total_h2 = total_nodes - total_h1 - total_byz
    total_honest = total_h1 + total_h2
    # honest_nodes_t1 = [0]*total_honest
    # honest_nodes_t2 = [0]*total_honest
    
    # for trial in range(trials):
    #     index = random.randint(0, total_honest-1)
    #     if index < total_h1:
    #         x = hyperormore_fast(total_nodes, total_h1 + total_byz, ksample, alpha)
    #         honest_nodes_t1[index] += x
    #         x = hyperormore_fast(total_nodes, total_h2, ksample, alpha)
    #         honest_nodes_t2[index] += x
    #     elif index >= total_h1:
    #         x = hyperormore_fast(total_nodes, total_h1, ksample, alpha)
    #         honest_nodes_t1[index] += x
    #         x = hyperormore_fast(total_nodes, total_h2 + total_byz, ksample, alpha)
    #         honest_nodes_t2[index] += x

    # print("Average H1 T1:", sum(honest_nodes_t1[:total_h1])/total_h1)
    # print("Average H1 T2:", sum(honest_nodes_t2[:total_h1])/total_h1)
    # print("Average H2 T1:", sum(honest_nodes_t1[total_h1:])/total_h2)
    # print("Average H2 T2:", sum(honest_nodes_t2[total_h1:])/total_h2)

    # # plt.plot(range(0, total_honest), honest_nodes_t1, '.', color='g')
    # # plt.plot(range(0, total_honest), honest_nodes_t2, '.')
    # # plt.show()
    h1t1 = trials * hyperormore_fast(total_nodes, total_h1 + total_byz, ksample, alpha)
    h1t2 = trials * hyperormore_fast(total_nodes, total_h2, ksample, alpha)
    h2t1 = trials * hyperormore_fast(total_nodes, total_h1, ksample, alpha)
    h2t2 = trials * hyperormore_fast(total_nodes, total_h2 + total_byz, ksample, alpha)
    data.append((h1t1/trials, h1t2/trials, h2t1/trials, h2t2/trials))
# plt.plot(ranges, [i[0] for i in data], label="h1t1")
# plt.plot(ranges, [i[1] for i in data], label="h1t2")

# plt.plot(ranges, [i[2] for i in data], label="h2t1")
# plt.plot(ranges, [i[3] for i in data], label="h2t2")
plt.semilogy(ranges, [i[0]/i[1] for i in data], label="h1t1/h1t2")
plt.semilogy(ranges, [i[1]/i[0] for i in data], label="h1t2/h1t1")
plt.semilogy(ranges, [i[3]/i[2] for i in data], label="h2t2/h2t1")
plt.semilogy(ranges, [i[2]/i[3] for i in data], label="h2t1/h2t2")
print("MIN h1t1/h1t2", min([i[0]/i[1] for i in data]))
print("MAX h1t1/h1t2", max([i[0]/i[1] for i in data]))
print("MIN h2t2/h2t1", min([i[3]/i[2] for i in data]))
print("MAX h2t2/h2t1", max([i[3]/i[2] for i in data]))
# plt.ylim((0, 100000000))
plt.legend()
plt.show()