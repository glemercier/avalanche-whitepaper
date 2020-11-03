import matplotlib.pyplot as plt

font=30
params = {'legend.fontsize': 40,
          'figure.figsize': (15, 5),
         'axes.labelsize': font,
         'axes.titlesize': font,
         'xtick.labelsize': font,
         'ytick.labelsize': font}
plt.rcParams.update(params)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Consecutive predicate used
# n = 100, b = 0, c = 100
# Total trials per node = 10**15
# Total trials = 10**15 * c
# k = 5, alpha = 4

# Consecutive = 200, to achieve red = 60, blue = 20 (maximal successes 80, since 60 + 20 byz = 80)
# phase shift = 17 (i.e. c/2 + 17 = 57)
# Security => 7.1e-9

# k = 10, alpha = 8
# Total trials per node = 10**15
# Total trials = 10**15 * c
# Consecutive = 150, to achieve red = 60, blue = 20 (maximal successes 80, since 60 + 20 byz = 80)
# phase shift = 13 (i.e. c/2 + 13 = 53)
# Security => 8.58e-35

# k = 20, alpha = 16
# Total trials per node = 10**15
# Total trials = 10**15 * c
# Consecutive = 120, to achieve red = 60, blue = 20 (maximal successes 80, since 60 + 20 byz = 80)
# phase shift = 11 (i.e. c/2 + 11 = 51)
# Security => 1.11e-99

# k = 40, alpha = 32
# Total trials per node = 10**15
# Total trials = 10**15 * c
# Consecutive = 110, to achieve red = 60, blue = 20 (maximal successes 80, since 60 + 20 byz = 80)
# phase shift = 11 (i.e. c/2 + 11 = 51)
# Security => 0

plt.subplot(2,1,1)
plt.semilogy([5, 10, 20, 40], [7.1e-9, 8.85e-35, 1.11e-99, 1.0e-300])
#plt.title("Varying k")
plt.ylabel(r"$p_f$")
plt.subplot(2,1,2)
plt.plot([5, 10, 20, 40], [200, 150, 120, 110])
plt.xlabel(r"$k$")
plt.ylabel(r"$l$")
plt.show()
