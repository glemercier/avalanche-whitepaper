import matplotlib.pyplot as plt
import pickle as pl
print("TRYING TO LOAD PICKLE")
fig = pl.load(open('output.pickle','rb'))
plt.show()