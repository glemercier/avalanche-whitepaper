import pykov
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.stats import hypergeom
import numpy as np

def hyper(N, S, k, alpha):
    rv = hypergeom(N, S, k)
    return 1 - rv.cdf(alpha) + rv.pmf(alpha)

def main():
    # N = 100
    # B = list(range(0, int(N/2), 1))
    # C1 = [int((N - i)/2) for i in B]
    # C2 = [N - B[i] - C1[i] for i in range(len(B))]
    # k = 20
    # alpha = list(range(int(k/2), k+1))
    N = 100
    B = 20
    C1 = 40
    C2 = 40
    k = 20
    alpha = 15
    states = []
    for i in range(int((C1+C2)/2), C1+C2+1):
        states.append((i, C1+C2-i))
    chain = pykov.Chain()    
    for i in range(1, len(states)-1):
        # print(states[i])
        # Prob of going up is prob choosing red and getting maj blue
        up =   (states[i][0]/(C1+C2))*hyper(N, states[i][1], k, alpha)
        # Prob of going down is prob of choosing blue
        down = (states[i][1]/(C1+C2))*hyper(N, states[i][0], k, alpha)
        # Prob of staying is prob of choosing blue and getting maj blue 
        # or choosing red and getting maj red
        # stay = ((states[i][0]/(C1+C2))*hyper(N, states[i][0], k, alpha))+((states[i][1]/(C1+C2))*hyper(N, states[i][1], k, alpha))
        stay = 1 - up - down
        chain[(str(states[i]),str(states[i-1]))] = up
        chain[(str(states[i]),str(states[i+1]))] = down
        chain[(str(states[i]),str(states[i]))] = stay
    chain[(str(states[0]),str(states[0]))] = 1 - (states[0][1]/(C1+C2))*hyper(N, states[0][0], k, alpha)
    chain[(str(states[0]),str(states[1]))] = (states[0][1]/(C1+C2))*hyper(N, states[0][0], k, alpha)
    chain[(str(states[len(states)-1]),str(states[len(states)-1]))] = 1 - (states[i][0]/(C1+C2))*hyper(N, states[i][1], k, alpha)
    chain[(str(states[len(states)-1]),str(states[len(states)-2]))] = (states[i][0]/(C1+C2))*hyper(N, states[i][1], k, alpha)
    # for state_i in states:
    #     for state_j in states:
    for k in chain:
        print(k, chain[k])
    p = pykov.Vector({str((40, 40)):1})
    print(chain.pow(p, 300000))
            

if __name__ == '__main__':
    main()
    