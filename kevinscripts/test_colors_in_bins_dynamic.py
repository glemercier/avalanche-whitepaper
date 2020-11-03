'''
In this test, we initialize a bin with RED and BLUE balls, initially zero. 
Then, for every step, we check the count, and if there's more RED balls than 
BLUE balls we add a BLUE ball with some probability p, we add a red with probability q, 
and add nothing with probability 1 - p - q (vice versa). It must be that p > q.
'''

import random
import matplotlib.pyplot as plt

def varying_rounds():
    for pp in [25]:
        averages = []
        distances = []
        rounds = range(50, 1100, 50)
        for r in rounds:
            trialsaverages = []
            trialsdistances = []
            p = pp
            q = 50 - p
            npq = 100 - p - q
            print(p, q, npq)
            for trial in range(1000):
                BINS = [0, 0]
                for _r in range(r):
                    if BINS[0] > BINS[1]:
                        index = random.randint(1, 100)
                        if 0 < index <= p:
                            BINS[1] += 1
                        elif p < index <= p + q:
                            BINS[0] += 1
                        else:
                            continue
                    elif BINS[1] > BINS[0]:
                        index = random.randint(1, 100)
                        if 0 < index <= p:
                            BINS[0] += 1
                        elif p < index <= p + q:
                            BINS[1] += 1
                        else:
                            continue
                    else:
                        index = random.randint(1, 100)
                        if index <= 50:
                            BINS[0] += 1
                        else:
                            BINS[1] += 1
                trialsaverages.append(abs(BINS[0]+BINS[1])/2)
                trialsdistances.append(abs(BINS[0]-BINS[1]))
            averages.append(sum(trialsaverages)/len(trialsaverages))
            distances.append(sum(trialsdistances)/len(trialsdistances))
        print(distances)
        plt.plot(rounds, distances, label="{}".format(p))
        # plt.plot(rounds, averages)
    plt.legend()
    plt.show()

def varying_p():
    averages = []
    distances = []
    pp = range(25, 51)
    for p in pp:
        trialsaverages = []
        trialsdistances = []
        q = 50 - p
        npq = 100 - p - q
        print(p, q, npq)
        for trial in range(50):
            BINS = [0, 0]
            for round in range(1000):
                if BINS[0] > BINS[1]:
                    index = random.randint(1, 100)
                    if 0 < index <= p:
                        BINS[1] += 1
                    elif p < index <= p + q:
                        BINS[0] += 1
                    else:
                        continue
                elif BINS[1] > BINS[0]:
                    index = random.randint(1, 100)
                    if 0 < index <= p:
                        BINS[0] += 1
                    elif p < index <= p + q:
                        BINS[1] += 1
                    else:
                        continue
                else:
                    index = random.randint(1, 100)
                    if index <= 50:
                        BINS[0] += 1
                    else:
                        BINS[1] += 1
            trialsaverages.append(abs(BINS[0]+BINS[1])/2)
            trialsdistances.append(abs(BINS[0]-BINS[1]))
        averages.append(sum(trialsaverages)/len(trialsaverages))
        distances.append(sum(trialsdistances)/len(trialsdistances))
    print(distances)
    plt.plot(pp, distances)
    # plt.plot(pp, averages)
    plt.show()

def main():
    # varying_p()
    varying_rounds()

if __name__ == '__main__':
    main()