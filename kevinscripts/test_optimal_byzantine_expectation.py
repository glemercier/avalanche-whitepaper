import random
import math
import itertools
from scipy.stats import hypergeom

def hyper(N, S, k, a):
    rv = hypergeom(N, S, k)
    return 1.0 - rv.cdf(a) + rv.pmf(a)

def main():
    _range = 10
    print("RANGE:", _range)
    outputs = []
    for seq in itertools.product("01", repeat=_range):
        correct = [-1]*3 + [1]*3
        byzantine = [-1]*3
        k = 5
        a = 3.75
        choices = [0]*len(correct)
        indeces = []
        # _range = int(3*len(correct)*math.log(len(correct), math.e))
        for i in range(_range):
            # Set byzantine
            if seq[i] == "0":
                byzantine = [-1]*3
            elif seq[i] == "1":
                byzantine = [1]*3
            # Select random correct node
            index = random.randint(0, len(correct)-1)
            indeces.append(index)
            choices[index] += 1
            correct_copy = correct[:]
            del correct_copy[index]
            sample = random.sample(correct_copy+byzantine, k)
            neg_count = len([i for i in sample if i < 0])
            pos_count = len([i for i in sample if i > 0])
            if neg_count >= a:
                # print(sample, neg_count, neg_count/len(sample))
                if correct[index] < 0:
                    correct[index] += -1
                elif correct[index] > 0:
                    correct[index] = -1
                else:
                    print("ERROR1")
                    exit(-1)
            elif pos_count >= a:
                # print(sample, pos_count, pos_count/len(sample))
                if correct[index] < 0:
                    correct[index] = 1
                elif correct[index] > 0:
                    correct[index] += 1
                else:
                    print("ERROR2")
                    exit(-1)
            elif pos_count < a and neg_count < a:
                continue
            else:
                print("ERROR3")
                exit(-1)

if __name__ == '__main__':
    main()