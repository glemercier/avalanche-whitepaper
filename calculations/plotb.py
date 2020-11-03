import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
from matplotlib import cm
from scipy.stats import hypergeom
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
import numpy as np

def hyper(N, S, k, alpha):
    rv = hypergeom(N, S, k)
    return 1 - rv.cdf(alpha) + rv.pmf(alpha)

def probOfStreak(numCoins, minHeads, headProb, saved=None):
    if saved == None: saved = {}
    ID = (numCoins, minHeads, headProb)
    if ID in saved: return saved[ID]
    else:
        if minHeads > numCoins or numCoins <= 0:
            result = 0
        else:
            result = headProb**minHeads
            for firstTail in range(1, minHeads+1):
                pr = probOfStreak(numCoins-firstTail, minHeads, headProb, saved)
                result += (headProb**(firstTail-1))*(1-headProb)*pr
        saved[ID] = result
        return result

def print_B_vs_alpha_liveness():
    N = 2000
    B = list(range(0, int(N/2), 1))
    C1 = [int((N - i)/2) for i in B]
    C2 = [N - B[i] - C1[i] for i in range(len(B))]
    k = 20
    alpha = list(range(int(k/2)+1, k+1))
    z = []
    for i in range(len(B)):
        print(i)
        zs = []
        for j in range(len(alpha)):
            r = hyper(N, N - B[i], k, alpha[j])
            # print("{},{},{}".format(B[i], C1[i], a), r)
            zs.append(r)
        z.append(zs)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y = np.meshgrid(alpha, B)
    z = np.array(z)
    print(x)
    print(y)
    print(z)
    ax.plot_surface(x, y, z, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    plt.ylabel("Byzantine")
    plt.xlabel("Alpha")
    plt.show()

def print_B_vs_alpha_safety():
    N = 500
    B = list(range(0, int(N/2), 1))
    C1 = [int((N - i)/2) for i in B]
    C2 = [N - B[i] - C1[i] for i in range(len(B))]
    k = 10
    alpha = list(range(int(k/2), k+1))
    zsafety = []
    zliveness = []
    cons = 20
    trials = 100
    for j in range(len(alpha)):
        # print(i)
        zssafety = []
        zsliveness = []
        for i in range(len(B)):    
            u = probOfStreak(trials, cons, hyper(N, B[i]+C1[i], k, alpha[j]))
            r_safety = (C1[i]**2)*(((1/(N-B[i]))*u)**2)
            v = 1 - probOfStreak(trials, cons, hyper(N, C1[i], k, alpha[j]))
            # r_safety = (len(C1)**2)*((probOfStreak(trials, cons, hyper(N, B[i]+C1[i], k, alpha[j])))**2)
            r_liveness = v
            # print("{},{},{}".format(B[i], C1[i], a), r)
            zssafety.append(r_safety)
            zsliveness.append(r_liveness)
            print(i, j, r_safety)                        
        zsafety.append(zssafety)
        zliveness.append(zsliveness)
    # fig = plt.figure()
    x, y = np.meshgrid(B, alpha)
    # max_z = max([max(i) for i in zsafety])
    # for i in range(len(zsafety)):
    #     for j in range(len(zsafety[i])):
    #         zsafety[i][j] /= max_z
    zsafety = np.array(zsafety)
    # max_z = max([max(i) for i in zliveness])
    # for i in range(len(zliveness)):
    #     for j in range(len(zliveness[i])):
    #         zliveness[i][j] /= max_z
    zliveness = np.array(zliveness)
    if True:
        fig, ax = plt.subplots(1, 1)
        pcm = ax.pcolor(x, y, zsafety, norm=colors.LogNorm(vmin=zsafety.min(), vmax=zsafety.max()),  cmap=cm.coolwarm)
        axim = ax.contour(x, y, zsafety, cmap=cm.Greys, norm = LogNorm())

        fmt = ticker.LogFormatterMathtext()
        fmt.create_dummy_axis()
        plt.clabel(axim, axim.levels, fmt=fmt)

        # CS = plt.contour(x, y, zsafety, norm=colors.LogNorm(vmin=zsafety.min(), vmax=zsafety.max()),  cmap=cm.binary)
        # fmt = ticker.LogFormatterMathtext()
        # fmt.create_dummy_axis()
        # plt.clabel(CS, CS.levels, fmt=fmt)
        # CS = plt.pcolor(x, y, zsafety, norm=colors.LogNorm(vmin=zsafety.min(), vmax=zsafety.max()),  cmap=cm.coolwarm)

        # plt.clabel(axim, axim.levels)
        fig.colorbar(pcm, ax=ax, extend='max')
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_surface(x, y, zsafety, cmap=cm.coolwarm,linewidth=0, antialiased=False)
        # ax = fig.add_subplot(212, projection='3d')
        # ax.plot_surface(x, y, zliveness,linewidth=0, antialiased=False)
        # ax.set_zlim(10e-10, 1)
    else:
        # Plot liveness
        fig, ax = plt.subplots(1, 1)
        pcm = ax.pcolor(x, y, zliveness, cmap=cm.coolwarm)
        axim = ax.contour(x, y, zliveness, cmap=cm.Greys)
        fmt = ticker.LogFormatterMathtext()
        fmt.create_dummy_axis()
        plt.clabel(axim, axim.levels, fmt=fmt)
        fig.colorbar(pcm, ax=ax, extend='max')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.ylabel(r"$\alpha k,\ k=10$")
    plt.xlabel("Byzantine |$\mathcal{B}$| size, ($|\mathcal{N}| = 500$)")
    plt.show()


def main():
    print_B_vs_alpha_safety()
    # print_B_vs_alpha_liveness()

if __name__ == '__main__':
    main()
    