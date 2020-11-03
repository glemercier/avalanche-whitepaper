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
    # T1globalpref = C1
    # T2globalpref = C2
    T1sum = mpm.mpf(0)
    T2sum = mpm.mpf(0)
    # gammas = []
    # xis = []
    gammas = mpm.mpf(0)
    xis = mpm.mpf(0)
    NUMTX = C*2

    C1nodes = [[mpm.mpf(0),mpm.mpf(0)] for _ in range(C1)]
    C2nodes = [[mpm.mpf(0),mpm.mpf(0)] for _ in range(C2)]
    
    history = HISTORY(5)

    log = []

    for i in range(0, NUMTX):
        # print(len(history.array))
        if i % 200 == 0:
            print(i, NUMTX)
        # print(i, T1sum, T2sum)
        # history.insert(0, [copy.deepcopy(C1nodes), copy.deepcopy(C2nodes)])
        # if len(history) == 6:
            # del history[len(history)-1]
        history.prepend(C1nodes, C2nodes)

        # _a = -len(history)
        # _b = len(history)
        # _sigma = 5

        # o = [trunc_normal(j, 0, _sigma, _a, _b) for j in range(_a, _b+1, 1)]
        # o = o[int(len(o)/2):]
        # for j in range(1, len(o)):
        #     o[j] = o[j] * 2

        T1globalpref = 0
        T2globalpref = 0
        for _i in range(len(history.array)):
            for __i in range(len(history.array[_i][0])):
                if history.array[_i][0][__i][0] > history.array[_i][0][__i][1]:
                    T1globalpref += (1*history.probabilities[_i])
                elif history.array[_i][0][__i][0] < history.array[_i][0][__i][1]:
                    T2globalpref += (1*history.probabilities[_i])
                else:
                    T1globalpref += (0.5*history.probabilities[_i])
                    T2globalpref += (0.5*history.probabilities[_i])
            for __i in range(len(history.array[_i][1])):
                if history.array[_i][1][__i][0] > history.array[_i][1][__i][1]:
                    T1globalpref += (1*history.probabilities[_i])
                elif history.array[_i][1][__i][0] < history.array[_i][1][__i][1]:
                    T2globalpref += (1*history.probabilities[_i])
                else:
                    T1globalpref += (0.5*history.probabilities[_i])
                    T2globalpref += (0.5*history.probabilities[_i])

        try:
            # assert mpm.mpf(C) >= T1sum + T2sum
            assert N >= (T1globalpref + B) - 1
            assert N >= (T2globalpref + B) - 1
        except AssertionError as e:
            print(T1sum, T2sum, T1globalpref, T2globalpref)
            raise e

        hypergamma1 = hyper(N, T1globalpref + B, k, alpha)
        hypergamma2 = hyper(N, T1globalpref, k, alpha)
        hyperxi1 = hyper(N, T2globalpref, k, alpha)
        hyperxi2 = hyper(N, T2globalpref + B, k, alpha)
        gamma1 = mpm.mpf(Z1/C)*hypergamma1
        gamma2 = mpm.mpf(Z2/C)*hypergamma2
        xi1 = mpm.mpf(Z1/C)*hyperxi1
        xi2 = mpm.mpf(Z2/C)*hyperxi2
        
        # gamma1 = mpm.mpf(Z1/C)*hyper(N, C1 + B + gammas - xis, k, alpha)
        # gamma2 = mpm.mpf(Z2/C)*hyper(N, C1 + gammas - xis, k, alpha)
        # xi1 = mpm.mpf(Z1/C)*hyper(N, C2 + xis - gammas, k, alpha)
        # xi2 = mpm.mpf(Z2/C)*hyper(N, C2 + B + xis - gammas, k, alpha)
        
        for j in range(len(C1nodes)):
            C1nodes[j][0] += (gamma1/C)
            C1nodes[j][1] += (xi1/C)
        for j in range(len(C2nodes)):
            C2nodes[j][0] += (gamma2/C)
            C2nodes[j][1] += (xi2/C)

        gammas += mpm.mpf(gamma2)
        xis += mpm.mpf(xi1)

        # print([(float(j[0]),float(j[1])) for j in C1nodes[:3]])
        # print([(float(j[0]),float(j[1])) for j in C2nodes[:3]])

        T1sum += (mpm.mpf(gamma1) + mpm.mpf(gamma2))
        T2sum += (mpm.mpf(xi1) + mpm.mpf(xi2))
        
        # print(id(C1nodes), [id(obj) for obj in [history.array[j][0] for j in range(len(history.array))]])
        # print(len(history.array), history.probabilities, history.length)
        if i % 1 == 0:
            logstr = "Time: {}; N, B, k, alpha = {}, {}, {}, {}; C1 = {}, C2 = {}\n".format(i, N, B, k, alpha, C1, C2)
            logstr += colored("T1globalpref="+str(T1globalpref)+";", "red") + " " + colored("T2globalpref="+str(T2globalpref)+";", "green") + "\n"
            logstr += colored("T1sum="+str(T1sum)+";", "red") + " " + colored("T2sum="+str(T2sum)+";", "green") + "\n"
            logstr += colored("C1.T1="+str(C1nodes[0][0])+";", "red") + " " + colored("C1.T2="+str(C1nodes[0][1])+";", "green") + " " + colored("C2.T1="+str(C2nodes[0][0])+";", "red") + " " + colored("C2.T2="+str(C2nodes[0][1])+";", "green") + "\n"
            # logstr += i, colored(float(gamma1), "red"), colored(float(gamma2), "red"), colored(float(xi1), "green"), colored(float(xi2), "green"))
            logstr += "----------------------------------\n"
            log.append(logstr)

        # T1globalpref += (mpm.mpf(gamma1) + mpm.mpf(gamma2))
        # T1globalpref -= (mpm.mpf(xi1) + mpm.mpf(xi2))
        # T2globalpref += (mpm.mpf(xi1) + mpm.mpf(xi2))        
        # T2globalpref -= (mpm.mpf(gamma1) + mpm.mpf(gamma2))
        # print(">> ", T1globalpref, T2globalpref)
        # print(":: ", T1sum, T2sum)
        # print(C1 + gammas - xis, C2 + xis - gammas)
        # print(T1sum, T1globalpref, ",", T2sum, T2globalpref)

    # print(gammas)
    # print(xis)
    T1globalpref = 0
    T2globalpref = 0
    for j in range(len(C1nodes)):
        if C1nodes[j][0] > C1nodes[j][1]:
            T1globalpref += 1
        elif C1nodes[j][0] < C1nodes[j][1]:
            T2globalpref += 1
        else:
            T1globalpref += 0.5
            T2globalpref += 0.5
    for j in range(len(C2nodes)):
        if C2nodes[j][0] > C2nodes[j][1]:
            T1globalpref += 1
        elif C2nodes[j][0] < C2nodes[j][1]:
            T2globalpref += 1
        else:
            T1globalpref += 0.5
            T2globalpref += 0.5
    print(t, ":::", T1sum, T1globalpref, ",", T2sum, T2globalpref)
    return (t[4], T1sum, T2sum, T1globalpref, T2globalpref, C1nodes[0], C2nodes[0], log)

def main():
    k = 20
    alpha = 14
    N = 500
    B = 100
    step = 25
    with Pool(8) as p:
        # r = p.map(exp, [(k, alpha, N, B, i) for i in range(1, int((N-B)/2)+step, step)])
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(121, 122, step)])
        # r = [exp((k, alpha, N, B, i)) for i in range(1, int((N-B)/2)+step, step)]
        # r = [exp((k, alpha, N, B, i)) for i in range(120, 122, step)]
    # r = [exp(200)]
    x = [i[0] for i in r]
    # y1 = [0 if math.isnan(i[1]) else i[1] for i in r]
    # y2 = [0 if math.isnan(i[2]) else i[2] for i in r]
    y1 = [i[1] for i in r]
    y2 = [i[2] for i in r]
    y3 = [i[3]/(N-B) for i in r]
    y4 = [i[4]/(N-B) for i in r]
    y5 = [i[5] for i in r]
    y6 = [i[6] for i in r]
    ylogs = [i[7] for i in r]
    print([len(i) for i in y5])
    fig = plt.figure(figsize=(15.0, 8.0))
    # plt.plot(x, y1, label="T1sum, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, y2, label="T2sum, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, y3, label="GT1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, y4, label="GT2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, [i[0] for i in y5], label="C1.T1")
    plt.plot(x, [i[1] for i in y5], label="C1.T2")
    plt.plot(x, [i[0] for i in y6], label="C2.T1")
    plt.plot(x, [i[1] for i in y6], label="C2.T2")
    plt.plot(x, [1/(i[1]/i[0]) for i in y5], label="C1.T2/C1.T1")
    plt.plot(x, [1/(i[1]/i[0]) for i in y6], label="C2.T2/C2.T1")
    # print("-------------------")
    # print(min(y1), max(y1))
    # print(min(y2), max(y2))
    # print(min(y3), max(y3))
    # print(min(y4), max(y4))
    with open("output.log", "w") as f:
        for log in ylogs:
            for entry in log:
                f.write(entry)
    # plt.legend()
    # plt.show()
    # plt.savefig("output.png")

if __name__ == '__main__':
    main()