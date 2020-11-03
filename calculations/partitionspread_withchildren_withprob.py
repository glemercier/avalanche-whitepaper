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

def binom(p, k, alpha):
    outcome = mpm.mpf(0)
    for i in range(alpha, k+1):
        outcome += ncr(k, i) * (p**i) * ((1 - p)**(k - i))
    return outcome

def exp(t):
    
    k = t[0]
    alpha = t[1]
    N = mpm.mpf(t[2])
    B = mpm.mpf(t[3])
    C = mpm.mpf(N - B)
    C1 = mpm.mpf(t[4])
    C2 = C - C1
    
    T1sum = mpm.mpf(0)
    T2sum = mpm.mpf(0)
    gammas = mpm.mpf(0)
    xis = mpm.mpf(0)
    
    NUMTX = int(C*2)

    # Each probs is P(T1>T2),   P(T2>T1),   P(T1==T2)
    C1Probs = [     mpm.mpf(1), mpm.mpf(0), mpm.mpf(0)]
    C2Probs = [     mpm.mpf(0), mpm.mpf(1), mpm.mpf(0)]

    # These are the expected sums in each bucket
    C1Exps = [mpm.mpf(0), mpm.mpf(0)]
    C2Exps = [mpm.mpf(0), mpm.mpf(0)]

    log = []

    T1globalpref = C1
    T2globalpref = C2

    T1globalprefhistory = []
    T2globalprefhistory = []

    CiProbshistory = []

    for i in range(0, 100):
        print(i, colored([float(_i) for _i in C1Probs], "red"), colored([float(_i) for _i in C2Probs], "green"))
            # if i % 200 == 0:
            #     print(i, NUMTX)
        
        CiProbshistory.append([C1Probs[:], C2Probs[:]])
        T1globalprefhistory.append(T1globalpref)
        T2globalprefhistory.append(T2globalpref)
        C1Probs_copy = copy.deepcopy(C1Probs)
        C2Probs_copy = copy.deepcopy(C2Probs)

        # For C1 bucket, probability that T1 > T2 after this round is equal 
        # to the probability that C1 gets selected, times the probability that in this 
        # round the bucket gets majority threshold T1.  
        # The probability that we get majority T1 for a single randomly chosen other bucket is 
        # probability that we choose a C1 bucket times its probability of majority + the
        # probability that we choose a C2 bucket times its probability of majority
        p = (C1/N)*C1Probs[0] + (C1/N)*C1Probs[2]*0.5 + (C2/N)*C2Probs[0] + (C2/N)*C2Probs[2]*0.5 + (B/N)*1
        alphaormoresuccesses = binom(p, k, alpha)
        C1Exps[0] += (1/C)*alphaormoresuccesses
        C1Probs_copy[0] = (C1/C)*alphaormoresuccesses

        # For C1 bucket, probability that T2 > T1 after this round is equal 
        # to the probability that C1 gets selected, times the probability that in this 
        # round the bucket gets majority threshold T2.
        # The probability that we get majority T2 for a single randomly chosen other bucket is 
        # probability that we choose a C1 bucket times its probability of majority + the
        # probability that we choose a C2 bucket times its probability of majority
        p = (C1/N)*C1Probs[1] + (C1/N)*C1Probs[2]*0.5 + (C2/N)*C2Probs[1] + (C2/N)*C2Probs[2]*0.5
        alphaormoresuccesses = binom(p, k, alpha)
        C1Exps[1] += (1/C)*alphaormoresuccesses
        C1Probs_copy[1] = (C1/C)*alphaormoresuccesses

        # Else
        C1Probs_copy[2] = 1 - C1Probs_copy[0] - C1Probs_copy[1]

        # --------------------------------------------------------------------------------

        # For C2 bucket, probability that T1 > T2 after this round is equal 
        # to the probability that C2 gets selected, times the probability that in this 
        # round the bucket gets majority threshold T1.  
        # The probability that we get majority T1 for a single randomly chosen other bucket is 
        # probability that we choose a C1 bucket times its probability of majority + the
        # probability that we choose a C2 bucket times its probability of majority
        p = (C1/N)*C1Probs[0] + (C1/N)*C1Probs[2]*0.5 + (C2/N)*C2Probs[0] + (C2/N)*C2Probs[2]*0.5
        alphaormoresuccesses = binom(p, k, alpha)
        C2Exps[0] += (1/C)*alphaormoresuccesses
        C2Probs_copy[0] = (C2/C)*alphaormoresuccesses

        # For C2 bucket, probability that T2 > T1 after this round is equal 
        # to the probability that C2 gets selected, times the probability that in this 
        # round the bucket gets majority threshold T2.
        # The probability that we get majority T2 for a single randomly chosen other bucket is 
        # probability that we choose a C1 bucket times its probability of majority + the
        # probability that we choose a C2 bucket times its probability of majority
        p = (C1/N)*C1Probs[1] + (C1/N)*C1Probs[2]*0.5 + (C2/N)*C2Probs[1] + (C2/N)*C2Probs[2]*0.5 + (B/N)*1
        print((C1/N)*C1Probs[1], (C1/N)*C1Probs[2]*0.5, (C2/N)*C2Probs[1], (C2/N)*C2Probs[2]*0.5, (B/N)*1)
        alphaormoresuccesses = binom(p, k, alpha)
        print(alphaormoresuccesses)
        C2Exps[1] += (1/C)*alphaormoresuccesses
        C2Probs_copy[1] = (C2/C)*alphaormoresuccesses

        # Else
        C2Probs_copy[2] = 1 - C2Probs_copy[0] - C2Probs_copy[1]

        C1Probs[0] = C1Probs_copy[0]
        C1Probs[1] = C1Probs_copy[1]
        C1Probs[2] = C1Probs_copy[2]
        C2Probs[0] = C2Probs_copy[0]
        C2Probs[1] = C2Probs_copy[1]
        C2Probs[2] = C2Probs_copy[2]

        # try:
        #     assert T1globalpref >= 0
        #     assert T2globalpref >= 0
        #     assert np.isclose(C, float(T1globalpref + T2globalpref))
        # except AssertionError as e:
        #     print("ASSERTION FAILED:", T1globalpref, T2globalpref, T1globalpref + T2globalpref)
        #     raise e

        # if i % 1 == 0:
        #     logstr = "Time: {}; N, B, k, alpha = {}, {}, {}, {}; C1 = {}, C2 = {}\n".format(i, N, B, k, alpha, C1, C2)
        #     logstr += colored("C1T1={};".format(gamma1), "red") + colored(" C2T1={};".format(gamma2), "red") + colored(" C1T2={};".format(xi1), "green") + colored(" C2T2={};\n".format(gamma2), "green")
        #     logstr += colored("T1globalpref="+str(T1globalpref)+";", "red") + " " + colored("T2globalpref="+str(T2globalpref)+";", "green") + "\n"
        #     logstr += colored("T1sum="+str(T1sum)+";", "red") + " " + colored("T2sum="+str(T2sum)+";", "green") + "\n"
        #     logstr += colored("C1.T1="+str(C1Probs[0])+";", "red") + " " + colored("C1.T2="+str(C1Probs[1])+";", "green") + " " + colored("C2.T1="+str(C2Probs[0])+";", "red") + " " + colored("C2.T2="+str(C2Probs[1])+";", "green") + "\n"
        #     # logstr += i, colored(float(gamma1), "red"), colored(float(gamma2), "red"), colored(float(xi1), "green"), colored(float(xi2), "green"))
        #     logstr += "----------------------------------\n"
        #     log.append(logstr)
        #     # print(logstr)

        # C1Probs[0] += (gamma1)
        # C1Probs[1] += (xi1)
        # C2Probs[0] += (gamma2)
        # C2Probs[1] += (xi2)
    
        # gammas += mpm.mpf(gamma2)
        # xis += mpm.mpf(xi1)

        # T1sum += (mpm.mpf(gamma1) + mpm.mpf(gamma2))
        # T2sum += (mpm.mpf(xi1) + mpm.mpf(xi2))
        
    # print(t, ":::", T1sum, T1globalpref, ",", T2sum, T2globalpref)
    return (t[4], C1Exps, C2Exps, T1globalpref, T2globalpref, C1Probs, C2Probs, log, T1globalprefhistory, T2globalprefhistory, CiProbshistory)

def main():
    fig = plt.figure(figsize=(15.0, 8.0))    
    k = 20
    alpha = 14
    N = 500
    B = 100
    step = 10
    with Pool(12) as p:
        # r = p.map(exp, [(k, alpha, N,  B, i) for i in range(1, int((N-B)), step)])
        # r = p.map(exp, [(k, alpha, N, B, i) for i in range(100, 300+step, step)])
        r = [exp((k, alpha, N, B, i)) for i in [200]]
        # r = [exp((k, alpha, N, B, i)) for i in range(100, 300+step, step)]
    x = [i[0] for i in r]
    y1 = [i[1] for i in r]
    y2 = [i[2] for i in r]
    y3 = [i[3]/(N-B) for i in r]
    y4 = [i[4]/(N-B) for i in r]
    y5 = [i[5] for i in r]
    y6 = [i[6] for i in r]
    ylogs = [i[7] for i in r]
    CiProbshistory = [i[10] for i in r][0]
    print(CiProbshistory)
    x = range(len(CiProbshistory))
    plt.plot(x, [i[0][0] for i in CiProbshistory], label="C1T1", marker="*")
    plt.plot(x, [i[0][1] for i in CiProbshistory], label="C1T2")
    plt.plot(x, [i[0][2] for i in CiProbshistory], label="C1None")
    plt.plot(x, [i[1][0] for i in CiProbshistory], label="C2T1")
    plt.plot(x, [i[1][1] for i in CiProbshistory], label="C2T2")
    plt.plot(x, [i[1][2] for i in CiProbshistory], label="C2None")
    # print([len(i) for i in y5])
    # plt.plot(x, [i[0] for i in y1], label="C1.T1 Exp, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, [i[1] for i in y1], label="C1.T2 Exp, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, [i[0] for i in y2], label="C2.T1 Exp, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, [i[1] for i in y2], label="C2.T2 Exp, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, y3, label="GT1, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, y4, label="GT2, k:{}, a:{}, N:{}, B:{}".format(k, alpha, N, B), linestyle='--')
    # plt.plot(x, [i[0] for i in y5], label="C1.T1.Prob")
    # plt.plot(x, [i[1] for i in y5], label="C1.T2.Prob")
    # plt.plot(x, [i[0] for i in y6], label="C2.T1.Prob")
    # plt.plot(x, [i[1] for i in y6], label="C2.T2.Prob")
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