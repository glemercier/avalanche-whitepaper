'''
Double absorbing states Snowflake (CTMC), 
mean-time to absorption to either 0 or N
using the matrix method (Tan method)
'''

from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import math
import mpmath as mpm
import numpy as np
# np.set_printoptions(threshold=np.nan)
from scipy.linalg import expm
from multiprocessing import Pool

mpm.mp.dps = 30

def hyper(pop, succ, sample, threshold):
    rv = hypergeom(pop, succ, sample)
    return mpm.mpf(1.0) - rv.cdf(threshold) + rv.pmf(threshold)

def print_matrix(matrix):
    print(np.array2string(matrix, max_line_width=np.inf))

def helper(tup):
    f, N, k, alpha, beta = tup
    # Now plot the prob of commit under bounded rounds (in this case, 
    # rounds is about c*log(N) where c is some constant)
    results = []
    print(f, N//2)
    correct = N - f

    births = []
    deaths = []

    '''
    State space = {alpha-1 blue,  alpha blue,  alpha+1 blue, ...,  N-alpha+1 blue}
                = {N-alpha+1 red, N-alpha red, N-alpha-1 red, ..., alpha-1 red}
    '''

    state_space = list(range(alpha-1, correct-alpha+2))

    for state in state_space:
        # Generate absorption states
        if state == alpha-1:
            births.append(0)
            deaths.append(0)
            continue
        # Generate transient states
        num_blue = state
        num_red = correct-state
        expected_blue_commit = num_blue * (hyper(N, num_blue, k, alpha)**beta)
        expected_red_commit = num_red * (hyper(N, num_red + f, k, alpha)**beta)
        if state < correct//2:
            # Birth rate (more blue) = Select red, sample blue majority
            sample_blue = hyper(N, num_blue+f, k, alpha)
            assert num_red >= expected_red_commit
            assert num_blue >= expected_blue_commit
            births.append((num_red-expected_red_commit)*sample_blue)
            # Death rate (less blue) = Select blue, sample red majority
            sample_red = hyper(N, num_red, k, alpha)
            deaths.append((num_blue-expected_blue_commit)*sample_red)
        elif state > correct//2:
            # Birth rate (more blue) = Select red, sample blue majority
            sample_blue = hyper(N, num_blue, k, alpha)
            assert num_red >= expected_red_commit
            assert num_blue >= expected_blue_commit
            births.append((num_red-expected_red_commit)*sample_blue)
            # Death rate (less blue) = Select blue, sample red majority
            sample_red = hyper(N, num_red+f, k, alpha)
            deaths.append((num_blue-expected_blue_commit)*sample_red)
        else:
            # Birth rate (more blue) = Select red, sample blue majority
            sample_blue = hyper(N, num_blue+f, k, alpha)
            assert num_red >= expected_red_commit
            assert num_blue >= expected_blue_commit
            births.append((num_red-expected_red_commit)*sample_blue)
            # Death rate (less blue) = Select blue, sample red majority
            sample_red = hyper(N, num_red+f, k, alpha)
            deaths.append((num_blue-expected_blue_commit)*sample_red)
    # print(state_space)
    normalized_state_space = list(range(0, len(state_space)))
    M = max(normalized_state_space)
    for start in range(M//2, M):
        print(start)
        numerator = mpm.mpf(0.0)
        for s in range(1, M):
            first_inner_sum = mpm.mpf(0.0)  
            for l in range(1, min(start, s) + 1):
                first_inner_product = mpm.mpf(1.0)
                for i in range(1, l):
                    first_inner_product *= deaths[i]
                second_inner_product = mpm.mpf(1.0)
                for j in range(l, s):
                    second_inner_product *= births[j]
                first_inner_sum += (first_inner_product * second_inner_product)
            second_inner_sum = mpm.mpf(0.0)
            for n in range(1, M - s - max(start - s, 0) + 1):
                first_inner_product = mpm.mpf(1.0)
                for i in range(s+1, M-n+1):
                    first_inner_product *= deaths[i]
                second_inner_product = mpm.mpf(1.0)
                for j in range(M-n+1, M):
                    second_inner_product *= births[j]
                second_inner_sum += (first_inner_product * second_inner_product)
            numerator += (first_inner_sum * second_inner_sum)
        denominator = mpm.mpf(0.0)
        for n in range(1, M+1):
            first_inner_product = mpm.mpf(1.0)
            for i in range(1, M-n+1):
                first_inner_product *= deaths[i]
            second_inner_product = mpm.mpf(1.0)
            for j in range(M-n+1, M):
                second_inner_product *= births[j]
            denominator += first_inner_product * second_inner_product
        results.append(numerator/denominator)
    return max(results)


PERCENT_NORMALIZED = True
MESSAGE_LOAD_NORMALIZED = False
def main():
    N = 200
    # for k, alpha, f in [(i, i, 0) for i in range(2, 8)]:
    byz_range = range(0, N//2, 1)
    with Pool(12) as p:
        results = p.map(helper, [(i, N, 10, 7, 150) for i in byz_range])
    #min_results = min(results)
    #results = [min_results/i for i in results]
    print([i / N for i in results])
    print([float(i) for i in results])
    plt.plot(byz_range, results)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
