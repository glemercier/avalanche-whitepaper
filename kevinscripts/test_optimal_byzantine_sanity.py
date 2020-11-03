import random
import math
import itertools
from scipy.stats import hypergeom

class CorrectNode(object):
    def __init__(self, color):
        self.color = color
        if self.color = 0:
            self.red_confidence = 1
            self.blue_confidence = 0
        else:
            self.blue_confidence = 1
            self.red_confidence = 0

    def get_color(self):
        return self.color                
        

def get_color(l):
    if abs(l[0]) > abs(l[1]):
        return -1
    elif abs(l[1] > abs(l[0])):
        return 1
    else:
        c = random.randint(0,1)
        if c:
            return 1
        else:
            return -1

def main():
    _range = 7
    outputs = []
    num_correct_each_color = 20
    num_byz = 15
    num_trials = 20
    k = 5
    a = 3.5
    outputs = []
    skipper = 0
    for byz_sequence in itertools.product("01", repeat=_range):
        for nature_sequence in itertools.product("01", repeat=_range):
            skipper += 1
            if skipper % 2 == 0:
                continue
            print("".join(byz_sequence), "".join(nature_sequence))
            trials = []
            for trial in range(num_trials):
                correct = [[-1,0]]*num_correct_each_color + [[0,1]]*num_correct_each_color
                byzantine = [-1]*num_byz            
                choices = [0]*len(correct)
                # Iterate over time
                for i in range(_range):
                    # print("--------------------------")
                    # print(byz_sequence)
                    # print(nature_sequence)
                    # print("TIME:", i, correct)
                    # Set byzantine
                    if byz_sequence[i] == "0":
                        byzantine = [[-1,0]]*num_byz
                    elif byz_sequence[i] == "1":
                        byzantine = [[0,1]]*num_byz
                    # Correct node color
                    corr_color = nature_sequence[i]
                    corr_color_index = None
                    if corr_color == "0":
                        # Red node was chosen, go ahead and find all red nodes 
                        # and select one of them uniformly at random
                        _indeces_red = []
                        for _i in range(len(correct)):
                            if get_color(correct[_i]) == -1:
                                _indeces_red.append(_i)
                        if _indeces_red:
                            # At least one red color node exists
                            corr_color_index = random.choice(_indeces_red)
                        else:
                            # print("No color red nodes exist, don't know what to do:")
                            # print(correct)
                            # print(_indeces_red)
                            # exit(-1)
                            continue
                    elif corr_color == "1":
                        # Blue node was chosen, go ahead and find all blue nodes 
                        # and select one of them uniformly at random
                        _indeces_blue = []
                        for _i in range(len(correct)):
                            if get_color(correct[_i]) == 1:
                                _indeces_blue.append(_i)
                        if _indeces_blue:
                            # At least one blue color node exists
                            corr_color_index = random.choice(_indeces_blue)
                        else:
                            # print("No color blue nodes exist, don't know what to do:")
                            # print(correct)                    
                            # exit(-1)
                            continue
                    correct_copy = correct[:]
                    del correct_copy[corr_color_index]
                    sample = random.sample(correct_copy+byzantine, k)
                    sample = [get_color(_i) for _i in sample]
                    neg_count = len([_i for _i in sample if _i == -1])
                    pos_count = len([_i for _i in sample if _i == 1])
                    if neg_count >= a:
                        correct[corr_color_index][0] += -1
                    elif pos_count >= a:
                        correct[corr_color_index][1] += 1
                    elif pos_count < a and neg_count < a:
                        continue
                pos_deltas = []
                neg_deltas = []
                for i in range(len(correct)):
                    if get_color(correct[i]) == -1:
                        neg_deltas.append(-1*(abs(correct[i][0])-abs(correct[i][1])))
                    if get_color(correct[i]) == 1:
                        pos_deltas.append(abs(correct[i][1]-correct[i][0]))
                if len(neg_deltas) == 0:
                    neg_deltas.append(0)
                if len(pos_deltas) == 0:
                    pos_deltas.append(0)
                trials.append([min(neg_deltas), max(pos_deltas)])
            distance = len([i for i in range(_range) if byz_sequence[i] != nature_sequence[i]])
            neg_av = sum(list(zip(*trials))[0])/len(trials)
            pos_av = sum(list(zip(*trials))[1])/len(trials)
            outputs.append([byz_sequence, nature_sequence, 
            distance, neg_av, pos_av, (neg_av+pos_av)/2, abs(neg_av-pos_av)])
    outputs = sorted(outputs, key=lambda x:(x[2]))
    
    import matplotlib.pyplot as plt    
    twoD = True
    if twoD:
        x = []
        z = []
        for i in outputs:
            print(i)
            x.append(i[2])
            z.append(i[5])
        plt.scatter(x, z)    
        plt.show()
    else:
        from mpl_toolkits.mplot3d import Axes3D    
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')    
        x = []
        y = []
        z = []
        for i in outputs:
            x.append(i[2])
            y.append(i[5])
            z.append(i[6])
        ax.scatter(x, y, z)
        plt.show()

if __name__ == '__main__':
    main()