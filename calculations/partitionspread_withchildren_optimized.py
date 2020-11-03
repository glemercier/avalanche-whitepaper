import random
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
import math
import mpmath as mpm
from multiprocessing import Pool
from scipy.special import comb
from termcolor import colored
import sys
import copy
import time
import numpy as np

mpm.mp.dps = 30

def ncr(n, r):
    return comb(n, r, exact=True)

def hyper(N, succ, sample, x):
    outcome = mpm.mpf(0)
    for i in range(x, sample+1):
        top = ncr(succ, i) * ncr(N - succ, sample - i)
        bottom = ncr(N, sample)
        outcome += (mpm.mpf(top)/mpm.mpf(bottom))
    return outcome

from scipy.special import erf

def trunc_normal(x, mu, sigma, a, b):
    return phi((x - mu)/sigma)/(sigma*(cdf((b-mu)/sigma) - cdf((a-mu)/sigma)))

def phi(eps):
    return (1/math.sqrt(2*math.pi))*math.pow(math.e, -0.5*math.pow(eps, 2))

def cdf(x):
    return 0.5*(1 + erf(x/math.sqrt(2)))

def normal(x, sigma):
    return (1/math.sqrt(2*math.pi*math.pow(sigma, 2)))*math.pow(math.e, -0.5*(x**2)*(1/math.pow(sigma, 2)))

def lineary(index, length):
    slope = -(2/math.pow(length, 2))
    return (slope*index) + (2/length)

class HISTORY():
    def __init__(self, length):
        self.array = []
        self.probabilities = None
        self.length = length
    def prepend(self, C1nodes, C2nodes):
        self.array.insert(0, [copy.deepcopy(C1nodes), copy.deepcopy(C2nodes)])
        if len(self.array) == self.length + 1:
            del self.array[self.length]
        self.probabilities = [lineary(i+0.5, len(self.array)) for i in range(len(self.array))]

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
    T1sum = mpm.mpf(0)
    T2sum = mpm.mpf(0)
    gammas = mpm.mpf(0)
    xis = mpm.mpf(0)
    NUMTX = C*2

    C1nodes = [mpm.mpf(0), mpm.mpf(0)]
    C2nodes = [mpm.mpf(0), mpm.mpf(0)]

    history = HISTORY(5)

    log = []

    T1globalpref = C1
    T2globalpref = C2

    T1globalprefhistory = []
    T2globalprefhistory = []

    for i in range(0, NUMTX):

        if i % 200 == 0:
            print(i, NUMTX)

        history.prepend(C1nodes, C2nodes)

        # T1globalpref = 0
        # T2globalpref = 0
        
        # for _i in range(len(history.array)):
        #     _C1T1 = history.array[_i][0][0]
        #     _C1T2 = history.array[_i][0][1]
        #     _C2T1 = history.array[_i][1][0]
        #     _C2T2 = history.array[_i][1][1]
        #     if _C1T1 > _C1T2:
        #         T1globalpref += (C1 * history.probabilities[_i])
        #     elif _C1T1 < _C1T2:
        #         T2globalpref += (C1 * history.probabilities[_i])
        #     else:
        #         T1globalpref += (0.5 * C1 * history.probabilities[_i])
        #         T2globalpref += (0.5 * C1 * history.probabilities[_i])
            
        #     if _C2T1 > _C2T2:
        #         T1globalpref += (C2 * history.probabilities[_i])
        #     elif _C2T1 < _C2T2:
        #         T2globalpref += (C2 * history.probabilities[_i])
        #     else:
        #         T1globalpref += (0.5 * C2 * history.probabilities[_i])
        #         T2globalpref += (0.5 * C2 * history.probabilities[_i])

        # T1globalpref = 0
        # T2globalpref = 0

        # _C1T1 = C1nodes[0]
        # _C1T2 = C1nodes[1]
        # _C2T1 = C2nodes[0]
        # _C2T2 = C2nodes[1]
        # if _C1T1 >= _C1T2:
        #     T1globalpref += (C1)*_C1T1
        # elif _C1T1 < _C1T2:
        #     T2globalpref += (C1)
        # # else:
        # #     T1globalpref += (0.5 * C1)
        # #     T2globalpref += (0.5 * C1)
        
        # if _C2T2 >= _C2T1:
        #     T2globalpref += (C2)
        # elif _C2T2 < _C2T1:
        #     T1globalpref += (C2)
        # # else:
        # #     T1globalpref += (0.5 * C2)
        # #     T2globalpref += (0.5 * C2)


        T1globalprefhistory.append(T1globalpref)
        T2globalprefhistory.append(T2globalpref)

        # Choose C1 and get T1
        gamma1 = mpm.mpf(C1/C)*hyper(N, T1globalpref + B, k, alpha)
        # Choose C2 and get T1
        gamma2 = mpm.mpf(C2/C)*hyper(N, T1globalpref, k, alpha)
        # Choose C1 and get T2
        xi1 = mpm.mpf(C1/C)*hyper(N, T2globalpref, k, alpha)
        # Choose C2 and get T2
        xi2 = mpm.mpf(C2/C)*hyper(N, T2globalpref + B, k, alpha)

        try:
            assert T1globalpref >= 0
            assert T2globalpref >= 0
            assert np.isclose(C, float(T1globalpref + T2globalpref))
        except AssertionError as e:
            print("ASSERTION FAILED:", T1globalpref, T2globalpref, T1globalpref + T2globalpref)
            raise e

        if i % 1 == 0:
            logstr = "Time: {}; N, B, k, alpha = {}, {}, {}, {}; C1 = {}, C2 = {}\n".format(i, N, B, k, alpha, C1, C2)
            logstr += colored("C1T1={};".format(gamma1), "red") + colored(" C2T1={};".format(gamma2), "red") + colored(" C1T2={};".format(xi1), "green") + colored(" C2T2={};\n".format(gamma2), "green")
            logstr += colored("T1globalpref="+str(T1globalpref)+";", "red") + " " + colored("T2globalpref="+str(T2globalpref)+";", "green") + "\n"
            logstr += colored("T1sum="+str(T1sum)+";", "red") + " " + colored("T2sum="+str(T2sum)+";", "green") + "\n"
            logstr += colored("C1.T1="+str(C1nodes[0])+";", "red") + " " + colored("C1.T2="+str(C1nodes[1])+";", "green") + " " + colored("C2.T1="+str(C2nodes[0])+";", "red") + " " + colored("C2.T2="+str(C2nodes[1])+";", "green") + "\n"
            # logstr += i, colored(float(gamma1), "red"), colored(float(gamma2), "red"), colored(float(xi1), "green"), colored(float(xi2), "green"))
            logstr += "----------------------------------\n"
            log.append(logstr)
            # print(logstr)

        C1nodes[0] += (gamma1)
        C1nodes[1] += (xi1)
        C2nodes[0] += (gamma2)
        C2nodes[1] += (xi2)
    
        gammas += mpm.mpf(gamma2)
        xis += mpm.mpf(xi1)

        T1sum += (mpm.mpf(gamma1) + mpm.mpf(gamma2))
        T2sum += (mpm.mpf(xi1) + mpm.mpf(xi2))
        
        # print(gamma1, xi1, gamma2, xi2)
        # T1globalpref += (-xi1 + gamma2)
        # T2globalpref += (-gamma2 + xi1)
        # T1globalpref_copy = T1globalpref
        # T2globalpref_copy = T2globalpref
        # T1globalpref_copy = T1globalpref + T1globalpref*(gamma1) - T1globalpref*(xi1) + T2globalpref*(gamma2) - T2globalpref*(xi2)
        T1globalpref_copy = T1globalpref - xi1 + gamma2
        # T2globalpref_copy = T2globalpref - T1globalpref*(gamma1) + T1globalpref*(xi1) - T2globalpref*(gamma2) + T2globalpref*(xi2)
        T2globalpref_copy = T2globalpref + xi1 - gamma2
        if T1globalpref_copy < 0:
            T1globalpref = 0
            T2globalpref = C
        elif T2globalpref_copy < 0:
            T1globalpref = C
            T2globalpref = 0
        else:
            T1globalpref = T1globalpref_copy
            T2globalpref = T2globalpref_copy

        
    # T1globalpref = 0
    # T2globalpref = 0
    # _C1T1 = C1nodes[0]
    # _C1T2 = C1nodes[1]
    # _C2T1 = C2nodes[0]
    # _C2T2 = C2nodes[1]
    # if _C1T1 > _C1T2:
    #     T1globalpref += (C1)
    # elif _C1T1 < _C1T2:
    #     T2globalpref += (C1)
    # else:
    #     T1globalpref += (0.5 * C1)
    #     T2globalpref += (0.5 * C1)
    
    # if _C2T1 > _C2T2:
    #     T1globalpref += (C2)
    # elif _C2T1 < _C2T2:
    #     T2globalpref += (C2)
    # else:
    #     T1globalpref += (0.5 * C2)
    #     T2globalpref += (0.5 * C2)
    print(t, ":::", T1sum, T1globalpref, ",", T2sum, T2globalpref)
    return (t[4], T1sum, T2sum, T1globalpref, T2globalpref, C1nodes, C2nodes, log, T1globalprefhistory, T2globalprefhistory)

def main():
    k = 20
    alpha = 14
    N = 500
    B = 100
    step = 10
    with Pool(12) as p:
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(1, int((N-B)), step)])
        # r = p.map(exp, [(k, alpha, N, B, i) fosr i in range(100, 300+step, step)])
        # r = [exp((k, alpha, N, B, i)) for i in [1]]
        # r = [exp((k, alpha, N, B, i)) for i in range(100, 300+step, step)]
    x = [i[0] for i in r]
    y1 = [i[1] for i in r]
    y2 = [i[2] for i in r]
    y3 = [i[3]/(N-B) for i in r]
    y4 = [i[4]/(N-B) for i in r]
    y5 = [i[5] for i in r]
    y6 = [i[6] for i in r]
    ylogs = [i[7] for i in r]
    # print([len(i) for i in y5])
    fig = plt.figure(figsize=(15.0, 8.0))
    # plt.plot(x, y1, label="T1sum, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, y2, label="T2sum, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, y3, label="GT1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, y4, label="GT2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, [i[0] for i in y5], label="C1.T1")
    plt.plot(x, [i[1] for i in y5], label="C1.T2")
    plt.plot(x, [i[0] for i in y6], label="C2.T1")
    plt.plot(x, [i[1] for i in y6], label="C2.T2")
    # print(y5)
    # print(y6)
    # plt.semilogy(x, [i[0]/i[1] for i in y5], label="C1.T1/C1.T2", linestyle="-.")
    # plt.semilogy(x, [i[0]/i[1] for i in y6], label="C2.T1/C2.T2")


    # # PLOTTING ON A SINGLE PARTITION, PER TIME BASIS
    # T1globalprefhistories = [i[8] for i in r]
    # T2globalprefhistories = [i[9] for i in r]
    # for i in range(len(T1globalprefhistories)):
    #     T1globalprefhistory = T1globalprefhistories[i]
    #     T2globalprefhistory = T2globalprefhistories[i]
    #     plt.plot(range(len(T1globalprefhistory)), T1globalprefhistory, label="T1")
    #     plt.plot(range(len(T2globalprefhistory)), T2globalprefhistory, label="T2")


    # with open("output_new.log", "w") as f:
    #     for log in ylogs:
    #         for entry in log:
    #             f.write(entry)
    plt.legend()
    plt.show()
    # plt.savefig("output_new.png")

if __name__ == '__main__':
    main()