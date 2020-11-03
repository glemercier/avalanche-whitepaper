import matplotlib.pyplot as plt
import numpy as np
from sympy import Symbol, Rational, binomial, expand_func
import mpmath as mp
import math
import sys
from multiprocessing import Pool
sys.setrecursionlimit(200000)
mp.mp.dps = 20000

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

def test1():
    n = 1000
    b = 200
    k = 10
    alpha = 8
    step = 100
    _range = range(step, n+step, step)
    trials = 400
    _successes = [int(i) for i in np.linspace(10, trials, num=5)]
    for successes in _successes:
        results = []
        for i in _range:
            print(i)
            h = hyper(n, i, k, alpha)
            results.append(binom(trials, successes, h))
        plt.semilogy(_range, results, label=successes)
    plt.plot(_range, [2**-32 for i in _range])
    plt.legend()
    plt.show()

def test2():
    _range = range(790, 810, 2)
    for trials in [140]:
        print(trials)
        results = []
        for i in _range:
            h = hyper(1000, i, 10, 8)
            results.append(binom(trials, int(trials*0.9), h))
        plt.semilogy(_range, results, label=trials)
    plt.plot(_range, [1/(10**math.log(2**(32), 10))]*len(_range))
    plt.plot(_range, [0.5]*len(_range))
    plt.legend()
    plt.show()

def test3():
    trials = 350
    results = []
    _range = [800]
    for i in _range:
        print(i)
        h = hyper(1000, i, 10, 8)
        results.append(binom(trials, int(0.90*trials), h))
    print("-------------------------------")
    print(results)
    for index, r in enumerate(results):
        if r >= math.pow(math.e, -44):
            print(_range[index], r)
    # plt.semilogy(_range, results)
    # plt.show()

def test4():
    plt.plot(range(0, 1000), [hyper(1000, i, 10, 8) for i in range(0, 1000)])
    plt.plot(range(0, 1000), [hyper(1000, i, 20, 16) for i in range(0, 1000)])
    plt.plot(range(0, 1000), [hyper(1000, i, 40, 32) for i in range(0, 1000)])
    plt.plot(range(0, 1000), [hyper(1000, i, 80, 64) for i in range(0, 1000)])
    plt.show()

def psi(lamb):
    return (2/mp.power(lamb, 2))*h(1 + lamb)

def h(lamb):
    return lamb*(mp.log(lamb, mp.e) - 1) + 1

def theorem3(n, D, N, lamb):
    n = mp.mp.mpf(n)
    D = mp.mp.mpf(D)
    N = mp.mp.mpf(N)
    lamb = mp.mp.mpf(lamb)
    mu = D/N
    sigma_square = mu*(1-mu)
    one_minus_fn = 1 - ((n - 1)/(N-1))
    p1 = mp.power(lamb, 2)/(2*sigma_square*one_minus_fn)
    p2 = lamb/(mp.sqrt(n)*sigma_square*one_minus_fn)
    return mp.power(mp.e, -1*p1*psi(p2))

def theorem3_modified(n, D, N):
    alpha = mp.mp.mpf(0.8)
    n = mp.mp.mpf(n)
    D = mp.mp.mpf(D)
    N = mp.mp.mpf(N)
    p = D/N
    assert p <= alpha
    lamb = alpha - p
    sigma_square = p*(1-p)
    one_minus_fn = 1 - ((n - 1)/(N-1))
    p1 = -1*mp.power(lamb, 2)/(2*sigma_square*one_minus_fn)
    p2 = lamb/(mp.sqrt(n)*sigma_square*one_minus_fn)
    return mp.sqrt(n)*mp.power(mp.e, p1*psi(p2))

def regularHoeffding(n, D, N):
    alpha = mp.mp.mpf(0.8)
    n = mp.mp.mpf(n)
    D = mp.mp.mpf(D)
    N = mp.mp.mpf(N)
    p = D/N
    t = alpha - p
    return mp.power(mp.power(mp.e, -2*mp.power(t, 2)*n), 1550)

def hoeffding(n, D, N, lamb):
    n = mp.mp.mpf(n)
    D = mp.mp.mpf(D)
    N = mp.mp.mpf(N)
    lamb = mp.mp.mpf(lamb)
    p = D/N
    t = (lamb/(n*mp.sqrt(n))) + (p/n) - p
    print(lamb, t)
    return mp.power(mp.e, -2*mp.power(t, 2)*n)

def test8():
    # Measure probability of at least 1 instance of T consecutive successes over P trials
    _range = list(range(200, 1000, 5))
    for T in [100]:
        for trials in [10**j for j in range(3, 20)]:
            results = []
            for i in _range:
                print(i)
                h = hyper(1000, i, 10, 8)
                results.append((trials-T)*mp.power(h, T))
            for index, r in enumerate(results):
                if r >= 0.9:
                    print(_range[index], r)
            plt.semilogy(_range, results, label=(trials, T))
    plt.plot(_range, [1/(10**math.log(2**(32), 10))]*len(_range))
    plt.plot(_range, [0.5]*len(_range))
    # plt.ylim((10**-100,1))
    plt.legend()
    plt.show()

def test7():
    # Measure probability of at least 1 instance of T consecutive successes over P trials
    _range = list(range(200, 1000, 5))
    for T in [10]:
        for trials in [1000]:
            results = []
            for i in _range:
                print(i)
                h = hyper(1000, i, 10, 8)
                phase1 = mp.mp.mpf(0.0)
                for j in range(0, T):
                    phase1 += (binomial(trials-T-1, j) * mp.power(h, j) * mp.power(mp.mp.mpf(1) - h, trials-T-1 - j))
                phase2 = mp.power(h, T)
                results.append(phase1*(1-h)*phase2)
            for index, r in enumerate(results):
                if r >= 0.9:
                    print(_range[index], r)
            plt.semilogy(_range, results, label=(trials, T))
    plt.plot(_range, [1/(10**math.log(2**(32), 10))]*len(_range))
    plt.plot(_range, [0.5]*len(_range))
    # plt.ylim((10**-100,1))
    plt.legend()
    plt.show()

def test5Varying():
    # Measure probability of at least 1 instance of the minimum size of T consecutive successes 
    # over P trials such that the probability is less than epsilon, over varying P
    n = 300
    c = 240
    b = 60
    epsilon = mp.power(2, -32)
    # _range = list(range(500, 1000, 1))
    results = []
    for trials in [10**i for i in range(1, 20)]:
        for i in range(trials):
            h = hyper(n, c, 10, 8)
            h = mp.power(h, i)
            if mp.mp.mpf(1) - mp.power(mp.mp.mpf(1)-h, trials-i+1) <= epsilon:
                results.append((trials, i))
                print(trials, i)
                break
    plt.semilogx(list(zip(*results))[0], list(zip(*results))[1])
    # plt.ylim((10**-100,1))
    plt.legend()
    plt.show()

def test5():
    # Measure probability of at least 1 instance of T consecutive successes over P trials
    n = 100
    c = 80
    b = 20
    _range = list(range(int(c/2), c, 1))
    # _range = list(range(500, 1000, 1))
    for T in [110]:
        for trials in [10**15]:
            results = []
            # with Pool(8) as pool:
            #     results = pool.map(prob_atleast_k_identical_cons_n_trials_p_probability, [(trials, T, i) for i in _range])
            for i in _range:
                print(i)
                h = hyper(n, i+b, 40, 32)
                h = mp.power(h, T)
                results.append(mp.mp.mpf(1) - mp.power(mp.mp.mpf(1)-h, trials-T+1))
            for index, r in enumerate(results):
                if r >= 0.9:
                    print(_range[index], r)
            plt.semilogy(_range, results, label=(trials, T))
    plt.plot(_range, [1/(10**math.log(2**(32), 10))]*len(_range))
    plt.plot(_range, [0.5]*len(_range))
    # plt.ylim((10**-100,1))
    plt.legend()
    plt.show()

def probOfStreak(numCoins, minHeads, headProb, saved=None):
	if saved == None: saved = {}
	ID = (numCoins, minHeads, headProb)

	if ID in saved:
		return saved[ID]
	else:
			if minHeads > numCoins or numCoins <= mp.mp.mpf(0):
				result = mp.mp.mpf(0)
			else:
				result = headProb ** minHeads
				for firstTail in range(1, minHeads+1):
					pr = probOfStreak(numCoins-firstTail, minHeads, headProb, saved)
					result += (headProb ** (firstTail-mp.mp.mpf(1))) * (mp.mp.mpf(1) - headProb) * pr
			saved[ID] = result
			return result

def prob_atleast_k_identical_cons_n_trials_p_probability(t):
    n, k, i = t
    print(n, k, i)
    p = hyper(1000, i, 10, 8)
    p = mp.mp.mpf(p)
    return probOfStreak(n, k, p)


def test6():
    # n = mp.mp.mpf(10)
    n = 20
    D = mp.mp.mpf(200)
    N = mp.mp.mpf(1000)
    # p = mp.mp.mpf(D/N)
    # ts = np.linspace(0.0, 0.5, 100)
    # lambs = np.linspace(0.0001, 30, 10000)
    # plt.semilogy(lambs, [theorem3(n, D, N, i) for i in lambs], label="Theorem 3")
    # plt.semilogy(lambs, [hoeffding(n, D, N, i) for i in lambs], label = "Hoeffding")
    # plt.semilogy(range(200, 800), [theorem3_modified(n, i, N) for i in (range(200, 800))], label = "Modified Theorem 3")
    _range = range(200, 1000)
    plt.semilogy(_range, [regularHoeffding(n, i, N) for i in _range], label = "Regular Hoeffding")
    plt.plot(_range, [1/(10**math.log(2**(64), 10))]*len(_range))
    plt.plot(_range, [0.5]*len(_range))
    plt.legend()
    plt.show()

def main():
    test1()
    # test2()
    # test3()
    # test4()
    # test5()
    # test6()
    # test7()
    # test8()
    # test5Varying()

if __name__ == '__main__':
    main()
