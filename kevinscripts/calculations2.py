import matplotlib.pyplot as plt
import numpy as np
from multiprocessing import Pool
from sympy import binomial, zeros, pprint
import mpmath as mp
import math
import time
mp.mp.dps = 20

def hyper(N, suc, sample, x):
    total = mp.mp.mpf(0)
    bottom = binomial(N, sample)
    for i in range(x, sample+1):
        top = binomial(suc, i) * binomial(N-suc, sample-i)
        total += (top/bottom)
    return total

def binom(trials, successes, p):
    total = mp.mp.mpf(0)
    p_ = mp.mp.mpf(p)
    for i in range(successes, trials+1):
        total += (binomial(trials, i) * mp.power(p_, i) * mp.power(1 - p_, trials - i))
    return total

def test3Helper(t):
    n, k, alpha, c, b = t
    for i in range(0, 205):
        # print(i, ((0.5*c + i)*hyper(n, 0.5*c + i - 1, k, alpha)), ((0.5*c-i)*hyper(n, 0.5*c - i + b - 1, k, alpha)))
        if (((0.5*c-i)/c)*hyper(n, 0.5*c + i, k, alpha)) >= (((0.5*c+i)/c)*hyper(n, 0.5*c - i + b, k, alpha)):
            return i

def test3( _range):
    r = []    
    with Pool(4) as p:
        r = p.map(test3Helper, _range)
    plt.plot([i[1] for i in _range], r, label="alpha = {}*k".format(_range[0][2]/_range[0][1]))

def test4Helper(t):
    n, k, alpha, c, b, consecutive = t
    foundI = None
    for i in range(0, 205):
        if ((0.5*c+i)*hyper(n, 0.5*c + i - 1, k, alpha)) >= ((0.5*c-i)*hyper(n, 0.5*c - i + b - 1, k, alpha)):
            print("Found i =", i, ";", n, k, alpha, c, b)
            foundI = i
            break
            # return mp.power(hyper(n, 0.5*c + i - 1, k, alpha), consecutive)
    for i in range(1, 100):
        if mp.power(hyper(n, 0.5*c + foundI - 1, k, alpha), i) <= mp.power(2, -32):
            return i*k

def test4( _range):
    r = []
    with Pool(4) as p:
        r = p.map(test4Helper, _range)
    plt.semilogy([i[1] for i in _range], r, label="alpha = {}*k".format(float(_range[0][2])/float(_range[0][1])))

def test1(n, k, alpha, c, b):
    # Plot the g(s_i), p(s_i), and q(s_i) of the MC produced by Theta_0
    r = []
    pis = []
    qis = []
    for i in range(int(c/2)):
        pi = ((0.5*c - i)/c)*hyper(n, 0.5*c + i, k, alpha)
        qi = ((0.5*c + i)/c)*hyper(n, 0.5*c - i + b, k, alpha)
        pis.append(pi)
        qis.append(qi)
        r.append(mp.power(pi/qi, -1))
    plt.plot([int(c/2)+i for i in range(int(c/2))], [1 for i in range(int(c/2))], color="black")
    plt.plot([int(c/2)+i for i in range(int(c/2))], [0.1 for i in range(int(c/2))], color="black")
    plt.plot([int(c/2)+i for i in range(int(c/2))], [0.01 for i in range(int(c/2))], color="black")
    plt.plot([int(c/2)+i for i in range(int(c/2))], [0.001 for i in range(int(c/2))], color="black")
    plt.plot([int(c/2)+i for i in range(int(c/2))], [0.0001 for i in range(int(c/2))], color="black")
    plt.plot([int(c/2)+i for i in range(int(c/2))], [0.00001 for i in range(int(c/2))], color="black")
    plt.plot([int(c/2)+i for i in range(int(c/2))], [0.000001 for i in range(int(c/2))], color="black")
    plt.semilogy([int(c/2)+i for i in range(int(c/2))], r, label="gis; k, a = {}, {}".format(k, alpha))
    plt.semilogy([int(c/2)+i for i in range(int(c/2))], pis, label="pis; k, a = {}, {}".format(k, alpha))
    plt.semilogy([int(c/2)+i for i in range(int(c/2))], qis, label="qis; k, a = {}, {}".format(k, alpha))

def birthDeathProbabilityNonEquivocation(n, k, alpha, c, b):
    P = zeros(c, c)
    for i in range(0, c-1):
        # Birth means node finds out about equivocation
        # Select node that doesn't know and have it know
        birth = ((c-i)/c)*hyper(n, i + 1, k, 1)
        # No deaths can exist
        P[i, i+1] = birth
        P[i, i] = mp.mp.mpf(1.0) - birth
    # Last guy finds out, stays there forever
    P[c-1, c-1] = mp.mp.mpf(1.0)
    for i in range(c):
        print(P.row(i), sum(P.row(i)))
    P = np.array(P).astype(np.float64)
    P = P**(100)
    # pprint(P.col(0))
    P = np.array(P).astype(np.float64)
    plt.matshow(P)
    plt.show()

def birthDeathProbabilityStartFromMinimalRequired(n, k, alpha, c, b):
    # First, find phase shift
    phaseShift = test3Helper((n, k, alpha, c, b))
    print(n, k, alpha, c, b)    
    print("Phase shift:", phaseShift)
    # Generate birth death probabilities
    P = zeros(int(c/2), int(c/2))
    for i in range(1, int(0.5*c)-1):
        # print(i)
        # Birth means a blue becomes red
        birth = ((0.5*c-i)/c)*hyper(n, 0.5*c + i, k, alpha)
        # Death means a red becomes blue
        death = ((0.5*c+i)/c)*hyper(n, 0.5*c - i + b, k, alpha)
        P[i, i-1] = death
        P[i, i+1] = birth
        P[i, i] = mp.mp.mpf(1.0) - birth - death
    # Half blue half red
    # Choose blue, get red
    P[0, 1] = (0.5)*hyper(n, c/2, k, alpha)
    # Nothing changes
    P[0, 0] = mp.mp.mpf(1.0) - (0.5)*hyper(n, c/2, k, alpha)
    # No blue, all red
    # Choose red, get blue
    P[int(0.5*c)-1, int(0.5*c)-2] = hyper(n, b, k, alpha)
    # Choose red, get nothing
    P[int(0.5*c)-1, int(0.5*c)-1] = mp.mp.mpf(1.0) - hyper(n, b, k, alpha)
    # for i in range(int(c/2)):
    #     print(P.row(i))
    # Truncate matrix with upper left being phase shift index
    Pprime = zeros(int(c/2)-b, int(c/2)-b)
    for i in range(0, int(c/2)-b):
        newRow = P.row(i+b)
        for j in range(0, int(c/2)-b):
            Pprime[i, j] = newRow[j+b]
    Pprime[0,0] = Pprime[0,0] + P[b,b-1]
    # print("----------------------------------------------")
    for i in range(int(c/2)-b):
        print(Pprime.row(i), sum(Pprime.row(i)))
    PprimeEigen = np.array(Pprime).astype(np.float64)
    w, v = np.linalg.eig(PprimeEigen)
    print(w)
    maxw = w[0]
    maxwindex = 0
    for i in range(1, len(w)):
        if w[i] > maxw:
            maxw = w[i]
            maxwindex = i
    print(maxwindex, v[maxwindex])
    print("===============================================")
    Pprime = Pprime**(10**15*c)
    # pprint(Pprime.col(0))
    Pprime = np.array(Pprime).astype(np.float64)
    plt.matshow(Pprime)
    plt.show()

def birthDeathProbabilityStartFromPhaseShift(n, k, alpha, c, b):
    # First, find phase shift
    phaseShift = test3Helper((n, k, alpha, c, b))
    print(n, k, alpha, c, b)    
    print("Phase shift:", phaseShift)
    # Generate birth death probabilities
    P = zeros(int(c/2), int(c/2))
    for i in range(1, int(0.5*c)-1):
        # print(i)
        # Birth means a blue becomes red
        birth = ((0.5*c-i)/c)*hyper(n, 0.5*c + i, k, alpha)
        # Death means a red becomes blue
        death = ((0.5*c+i)/c)*hyper(n, 0.5*c - i + b, k, alpha)
        P[i, i-1] = death
        P[i, i+1] = birth
        P[i, i] = mp.mp.mpf(1.0) - birth - death
    # Half blue half red
    # Choose blue, get red
    P[0, 1] = (0.5)*hyper(n, c/2, k, alpha)
    # Nothing changes
    P[0, 0] = mp.mp.mpf(1.0) - (0.5)*hyper(n, c/2, k, alpha)
    # No blue, all red
    # Choose red, get blue
    P[int(0.5*c)-1, int(0.5*c)-2] = hyper(n, b, k, alpha)
    # Choose red, get nothing
    P[int(0.5*c)-1, int(0.5*c)-1] = mp.mp.mpf(1.0) - hyper(n, b, k, alpha)
    # for i in range(int(c/2)):
    #     print(P.row(i))
    # Truncate matrix with upper left being phase shift index
    Pprime = zeros(int(c/2)-phaseShift, int(c/2)-phaseShift)
    for i in range(0, int(c/2)-phaseShift):
        newRow = P.row(i+phaseShift)
        for j in range(0, int(c/2)-phaseShift):
            Pprime[i, j] = newRow[j+phaseShift]
    Pprime[0,0] = Pprime[0,0] + P[phaseShift,phaseShift-1]
    # print("----------------------------------------------")
    for i in range(int(c/2)-phaseShift):
        print(Pprime.row(i), sum(Pprime.row(i)))
    PprimeEigen = np.array(Pprime).astype(np.float64)
    w, v = np.linalg.eig(PprimeEigen)
    print(w)
    maxw = w[0]
    maxwindex = 0
    for i in range(1, len(w)):
        if w[i] > maxw:
            maxw = w[i]
            maxwindex = i
    print(maxwindex, v[maxwindex])
    print("===============================================")
    Pprime = Pprime**100000
    # pprint(Pprime.col(0))
    Pprime = np.array(Pprime).astype(np.float64)
    plt.matshow(Pprime)
    plt.show()

def Bn(pi, qi, n):
    bn = mp.mp.mpf(0.0)
    for x in range(n):
        piprod = mp.mp.mpf(1.0)
        qiprod = mp.mp.mpf(1.0)
        if x == 0:
            bn += mp.mp.mpf(1.0)
        else:
            for j in range(0, x):
                piprod = piprod * pi[j]
            for j in range(1, x+1):
                qiprod = qiprod * qi[j]
            bn += (piprod/qiprod)
    return bn

def birthDeathChainInvariantPDF(n, k, alpha, c, b):
    P = zeros(c+1, c+1)
    for i in range(1, c):
        print(i)
        # Birth means a blue becomes red
        birth = ((c-i)/c)*hyper(n, i, k, alpha)
        # Death means a red becomes blue
        death = (i/c)*hyper(n, c - i + b, k, alpha)
        P[i, i-1] = death
        P[i, i+1] = birth
        P[i, i] = mp.mp.mpf(1.0) - birth - death
    # Choose blue, get red
    P[0, 1] = hyper(n, 0, k, alpha)
    # Nothing changes
    P[0, 0] = mp.mp.mpf(1.0) - hyper(n, 0, k, alpha)
    # No blue, all red
    # Choose red, get blue
    P[c, c-1] = hyper(n, b, k, alpha)
    # Choose red, get nothing
    P[c, c] = mp.mp.mpf(1.0) - hyper(n, b, k, alpha)

    pi = []
    qi = []
    pi.append(P[0, 1])
    qi.append(mp.mp.mpf(0.0))
    for i in range(1, c):
        pi.append(P[i, i+1])
        qi.append(P[i, i-1])
    pi.append(mp.mp.mpf(0.0))
    qi.append(P[c, c-1])

    invariant = []
    for x in range(2, c+1):
        piprod = mp.mp.mpf(1.0)
        qiprod = mp.mp.mpf(1.0)
        for i in range(0, x):
            piprod = piprod*pi[i]
        for i in range(1, x+1):
            qiprod = qiprod*qi[i]
        print(piprod, qiprod)
        invariant.append((1/Bn(pi, qi, c))*(piprod/qiprod))
    plt.plot(range(2, c+1), invariant)
    plt.show()

def birthDeathProbabilityEntireRange(n, k, alpha, c, b):
    phaseShift = test3Helper((n, k, alpha, c, b))
    print(n, k, alpha, c, b)    
    print("Phase shift:", phaseShift)
    P = zeros(c, c)
    for i in range(1, c-1):
        print(i)
        # Birth means a blue becomes red
        birth = ((c-i)/c)*hyper(n, i, k, alpha)
        # Death means a red becomes blue
        death = (i/c)*hyper(n, c - i + b, k, alpha)
        P[i, i-1] = death
        P[i, i+1] = birth
        P[i, i] = mp.mp.mpf(1.0) - birth - death
    # Choose blue, get red
    P[0, 1] = hyper(n, 0, k, alpha)
    # Nothing changes
    P[0, 0] = mp.mp.mpf(1.0) - hyper(n, 0, k, alpha)
    # No blue, all red
    # Choose red, get blue
    P[c-1, c-2] = hyper(n, b, k, alpha)
    # Choose red, get nothing
    P[c-1, c-1] = mp.mp.mpf(1.0) - hyper(n, b, k, alpha)
    for i in range(c):
        print(P.row(i), sum(P.row(i)))
    Pprime = np.array(P).astype(np.float64)
    w, v = np.linalg.eig(Pprime)
    print(w)
    maxw = w[0]
    maxwindex = 0
    for i in range(1, len(w)):
        if w[i] > maxw:
            maxw = w[i]
            maxwindex = i
    print(maxwindex, v[maxwindex])
    for i in range(10, 100, 10):
        print((P**i)[60, 57])
    # for i in range(6):
    #     print(P[70, 57])
    #     P = P**10
    P = np.array(P).astype(np.float64)
    plt.matshow(P)
    plt.show()


def birthDeathProbability(n, k, alpha, c, b):
    # bi = []
    # for i in range(0, c-1):
    #     bi.append(((c-i)/c)*hyper(n, i - 1, k, alpha))
    # bi.append(mp.mp.mpf(0.0))
    # di = [mp.mp.mpf(0.0)]
    # for i in range(1, c):
    #     di.append(((i)/c)*hyper(n, c - i + b - 1, k, alpha))
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
    P = P
    for i in range(int(c/2)):
        print(i, P.col(i))
    print("=======================================================")
    print("=======================================================")
    P = np.array(P).astype(np.float64)
    plt.contour(P)
    plt.legend()
    plt.show()
    # pprint(P**10000)

if __name__=="__main__":
    if True:
        i = 20
        # birthDeathProbability(10*i, 20, 16, 8*i, 2*i)
        # birthDeathProbabilityStartFromPhaseShift(10*i, 5, 4, 8*i, 2*i)
        # birthDeathProbabilityStartFromMinimalRequired(10*i, 40, 32, 8*i, 2*i)
        birthDeathProbabilityNonEquivocation(10*i, 10, 8, 8*i, 2*i)
        # birthDeathProbabilityEntireRange(10*i, 5, 4, 8*i, 2*i)
        # birthDeathChainInvariantPDF(10*i, 10, 8, 8*i, 2*i)
        exit()
    elif False:
        i = 20
        # test1(1000, 1, 1, 800, 200)
        # test1(1000, 5, 4, 800, 200)
        # test1(1000, 10, 8, 800, 200)
        test1(10*i, 20, 16, 8*i, 2*i)
        # test1(1000, 40, 32, 800, 200)
    elif True:
        step = 1
        _range = [(1000, i, i, 800, 200) for i in range(1, 100+step, step)]
        test3(_range)
    elif False:
        step = 1
        # _range = [(1000, i, int(0.8*i), 800, 200, 1) for i in range(1, 40+step, step)]
        _range = [(1000, i, i, 800, 200, 1) for i in range(1, 40+step, step)]
        test4( _range)
    plt.legend()
    plt.show()