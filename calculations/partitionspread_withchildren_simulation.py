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
        self.probabilities = []
        self.probabilities_normalized = []
        self.length = length
        
    def prepend(self, nodes):
        self.array.insert(0, copy.deepcopy(nodes))
        if len(self.array) == self.length + 1:
            del self.array[self.length]
        if len(self.probabilities) < self.length:
            self.probabilities = [lineary(i+0.5, len(self.array)) for i in range(len(self.array))]
            self.probabilities_normalized = [self.probabilities[0]]
            for i in range(1, len(self.probabilities)):
                self.probabilities_normalized.append(self.probabilities[i] + self.probabilities_normalized[i-1])
            self.probabilities_normalized = [int(i*10000) for i in self.probabilities_normalized]


    def select_random_history(self):
        # print("Probabilities:", self.probabilities)
        index = random.randint(0, self.probabilities_normalized[len(self.probabilities_normalized)-1])
        for i in range(len(self.probabilities_normalized)):
            if self.probabilities_normalized[i] >= index:
                return i
        return -1
        

class Node():
    def __init__(self, index, identity, byzantine, k, alpha):
        self.index = index
        self.identity = identity
        self.transactions = {identity:0}
        self.byzantine = byzantine
        self.sample = []
        self.indeces_sample = set()
        self.k = k
        self.alpha = alpha

    def pretty_print(self):
        if "T1" not in self.transactions:
            return "_,"+str(self.transactions["T2"])
        if "T2" not in self.transactions:
            return str(self.transactions["T1"])+",_"
        else:
            return str(self.transactions["T1"])+","+str(self.transactions["T2"])
    
    def get_confidences(self):
        if "T1" not in self.transactions:
            return (0, self.transactions["T2"])
        if "T2" not in self.transactions:
            return (self.transactions["T1"], 0)
        else:
            return (self.transactions["T1"], self.transactions["T2"])

    def onpoll(self, identity):
        # Byzantine node
        if self.byzantine:
            return identity

        # Else, honest node
        # If this node doesn't even know about T1, just return T2
        if "T1" not in self.transactions:
            return "T2"
        # If this node doesn't even know about T2, just return T1
        elif "T2" not in self.transactions:
            return "T1"
        
        # Both transactions available
        # If T1 is bigger than T2
        if self.transactions["T1"] > self.transactions["T2"]:
            return "T1"
        # If T2 is more preferred
        elif self.transactions["T2"] > self.transactions["T1"]:
            return "T2"
    
        # Transactions equal
        if random.randint(0,1):
            return "T2"
        return "T1"
        
    # Polls a single node from the history
    def poll(self, history, nodes):
        # print(self.index, self.identity, len(self.sample))
        # First, select some time in history:
        hindex = history.select_random_history()
        if hindex == -1:
            print("Something went wrong...")
            exit(-1)
        # Now, select random node from history
        nindex = -1
        while True:
            nindex = random.randint(0, len(nodes)-1)
            if nindex == self.index:
                continue
            if nindex in self.indeces_sample:
                continue
            # Good index
            break
        # print("History index: {}, node index: {}".format(hindex, nindex))
        # print(len(history.array), len(history.array[hindex]))
        preferred_transaction = history.array[hindex][nindex].onpoll(self.identity)
        self.sample.append(preferred_transaction)
        self.indeces_sample.add(nindex)
        # If sample is now full, update:
        if len(self.sample) == self.k:
            t1count = self.sample.count("T1")
            t2count = self.sample.count("T2")
            assert self.k == t1count + t2count
            if (t1count > 0) and ("T1" not in self.transactions):
                self.transactions["T1"] = 0
            elif (t2count > 0) and ("T2" not in self.transactions):
                self.transactions["T2"] = 0
            if t1count >= self.alpha:
                self.transactions["T1"] += 1
            elif t2count >= self.alpha:
                self.transactions["T2"] += 1
            # Reset sample now
            self.sample = []
            self.indeces_sample = set()

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
    nodes = []
    for i in range(C1):
        nodes.append(Node(i, "T1", False, k, alpha))
    for i in range(C2):
        nodes.append(Node(i+C1, "T2", False, k, alpha))
    for i in range(B):
        nodes.append(Node(0, "NONE", True, k, alpha))

    history = HISTORY(100)

    log = []

    for i in range(0, 1000):
        # if i % 200 == 0:
        print(i)
        history.prepend(nodes)
        # Select random node to execute poll
        for n in range(random.randint(100, 1000)):
            random_node = random.randint(0, C1+C2-1)
            nodes[random_node].poll(history, nodes)
    
    C1T1confidences = 0
    C1T2confidences = 0
    C2T1confidences = 0
    C2T2confidences = 0
    # Dump to file: 
    # with open("nodes_k{}_alpha{}_B{}_C1{}_C2{}".format(k, alpha, B, C1, C2), "w+") as f:
    st = "C1 Nodes ----------------------------\n"
    for i in range(C1):
        st += "{}|".format(nodes[i].pretty_print())
        C1T1confidences += nodes[i].get_confidences()[0]
        C1T2confidences += nodes[i].get_confidences()[1]
    st += "\nC2 Nodes ----------------------------\n"
    for i in range(C1, C1+C2):
        st += "{}|".format(nodes[i].pretty_print())
        C2T1confidences += nodes[i].get_confidences()[0]
        C2T2confidences += nodes[i].get_confidences()[1]
    print(st)
    return (C1, C1T1confidences/C1, C1T2confidences/C1, C2T1confidences/C2, C2T2confidences/C2)

def main():
    k = 20
    alpha = 14
    N = 500
    B = 100
    step = 10
    with Pool(12) as p:
        # r = p.map(exp, [(k, alpha, N, B, i) for i in range(1, int((N-B))+step, step)])
        r = p.map(exp, [(k, alpha, N, B, i) for i in range(100, 300+step, step)])
        # r = [exp((k, alpha, N, B, i)) for i in range(1, int((N-B)/2)+step, step)]
        # r = [exp((k, alpha, N, B, i)) for i in range(120, 121, step)]
    x = [i[0] for i in r]
    y1 = [i[1] for i in r]
    y2 = [i[2] for i in r]
    y3 = [i[3] for i in r]
    y4 = [i[4] for i in r]
    fig = plt.figure(figsize=(15.0, 8.0))
    plt.plot(x, y1, label="C1T1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, y2, label="C1T2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, y3, label="C2T1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    plt.plot(x, y4, label="C2T2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, [i[0]/i[1] for i in y5], label="C1.T2/C1.T1")
    # plt.plot(x, [i[0]/i[1] for i in y6], label="C2.T2/C2.T1")
    # with open("output.log", "w") as f:
    #     for log in ylogs:
    #         for entry in log:
    #             f.write(entry)
    plt.legend()
    plt.show()
    # plt.savefig("output.png")

if __name__ == '__main__':
    main()