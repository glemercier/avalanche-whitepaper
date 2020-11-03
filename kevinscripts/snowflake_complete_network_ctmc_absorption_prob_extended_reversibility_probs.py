from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import math
import mpmath as mpm
from multiprocessing import Pool

mpm.mp.dps = 1000

def hyper(pop, succ, sample, threshold):
    rv = hypergeom(pop, succ, sample)
    return mpm.mpf(1.0) - rv.cdf(threshold) + rv.pmf(threshold)

def helper(_input):
    f, N, k, alpha, beta = _input
    correct = N - f
    births = []
    deaths = []

    '''
    State space = {alpha-1 blue,  alpha blue,  alpha+1 blue, ...,  N-alpha+1 blue}
                = {N-alpha+1 red, N-alpha red, N-alpha-1 red, ..., alpha-1 red}
    '''

    state_space = list(range(alpha-1, correct+1))

    for state in state_space:
        # Generate absorption states
        if state == alpha-1:
            births.append(0)
            deaths.append(0)
            continue
        if state == correct:
            births.append(0)
            deaths.append(0)
            continue
        # Generate transient states
        num_blue = state
        num_red = correct-state
        expected_blue_commit = num_blue * (hyper(N, num_blue, k, alpha)**beta)
        expected_red_commit = num_red * (hyper(N, num_red + f, k, alpha)**beta)
        # print(expected_blue_commit, expected_red_commit)
        # Birth rate (more blue) = Select red, sample blue majority
        sample_blue = hyper(N, num_blue, k, alpha)
        births.append((num_red-expected_red_commit)*sample_blue)
        # Death rate (less blue) = Select blue, sample red majority
        sample_red = hyper(N, num_red + f, k, alpha)
        deaths.append((num_blue-expected_blue_commit)*sample_red)
        # deaths.append((num_blue-1)*sample_red)

    # print(state_space)
    normalized_state_space = list(range(0, len(state_space)))
    # print(normalized_state_space)
    # print(births)
    # print(deaths)
    # exit()

    M = max(normalized_state_space)
    print(M)
    results_probs = []
    for s in range(M):
        print(s)
        # Calculate absorption into all-red state from state s

        # Calculate numerator
        numerator_total = mpm.mpf(0.0)

        for l in range(1, M-s+1):
            # print(l, list(range(1, M-l+1)), list(range(M-l+1, M-1+1)))
            prod1 = mpm.mpf(1.0)
            for i in list(range(1, M-l+1)):
                prod1 *= deaths[i]
            prod2 = mpm.mpf(1.0)
            for j in list(range(M-l+1, M-1+1)):
                prod2 *= births[j]
            numerator_total += prod1*prod2

        # Calculate denominator
        denominator_total = mpm.mpf(0.0)
        for l in range(1, M+1):
            # print(l, list(range(1, M-l+1)), list(range(M-l+1, M-1+1)))
            prod1 = mpm.mpf(1.0)
            for i in list(range(1, M-l+1)):
                prod1 *= deaths[i]
            prod2 = mpm.mpf(1.0)
            for j in list(range(M-l+1, M-1+1)):
                prod2 *= births[j]
            denominator_total += prod1*prod2

        p = hyper(N, s + alpha + f, k, alpha)**beta
        # print(N, s + alpha + f, k, alpha, beta, p)
        prob_reverse = numerator_total/denominator_total

        # # Prob single-commit
        # prob_single_commit_at_starting_point_s = correct*p*((mpm.mpf(1.0)-p)**(correct-1))
        # print("::", prob_reverse, prob_single_commit_at_starting_point_s)
        # prob_safety_failure = prob_single_commit_at_starting_point_s * prob_reverse

        # Prob at least single-commit
        prob_at_least_single_commit_at_starting_point_s = mpm.mpf(1.0) - (mpm.mpf(1.0) - p)**correct
        # print("::", prob_reverse, prob_at_least_single_commit_at_starting_point_s)
        prob_safety_failure = prob_at_least_single_commit_at_starting_point_s * prob_reverse

        # print(N, s, f, prob_safety_failure)
        results_probs.append(prob_safety_failure)
    # plt.semilogy(state_space[:-1], results_probs)
    # plt.show()
    # plt.plot(state_space[:-1], results_probs)
    # plt.show()
    return max(results_probs)

def experiment(N, k, alpha, beta):
    byz = range(0, N//2, 1)
    results = []
    # helper((40, 200, 30, 30, beta))
    # exit()
    with Pool(12) as p:
        results = p.map(helper, [(f, N, k, alpha, beta) for f in byz])
    print([float(i) for i in results])
    plt.plot([int(100*(f/N)) for f in byz], results)
    # plt.plot(byz, results)

def main():
    N = 200
    for k, alpha, beta in [(10, 7, 150)]:
        experiment(N, k, alpha, beta)
    plt.legend()
    plt.show()

if __name__ == '__main__':
    main()
