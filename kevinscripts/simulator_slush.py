import random
import numpy as np

def slush(sizen):
    k = 10
    alpha = 8
    iters = 100
    trials = 0
    numbers = []
    for trial in range(iters):
        n = int(sizen/2)*[0] + int(sizen/2)*[1]    
        i = 1
        while True:
            rcount = n.count(0)
            if rcount == sizen or rcount == 0:
                # print(i)
                trials +=i
                numbers.append(i/sizen)
                break
            u = random.randint(0, sizen-1)
            sample = [n[j] for j in random.sample(range(0, sizen), k)]
            if sample.count(0) >= alpha:
                n[u] = 0
            elif sample.count(1) >= alpha:
                n[u] = 1
            i += 1
    print(sizen, (trials/iters)/sizen)
    print(np.std(numbers))
    print('-------------------------')

def main():
    for nsize in [600, 1200, 2400, 4800, 9600]:
        slush(nsize)

if __name__ == '__main__':
    main()