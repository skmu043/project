import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)
shel = shelve.open(data_dr + "/" + str(data_archives[0]) + "/dyke.truncation.rk4.data")

ENV_START_ = []
ALIVE_THREASHOLD = 0

number_alive_diff = []
RANGE_R = 0

alive_threshholds=[]

for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.truncation.rk4.data")

    try:

        number_alive_global_start = s['number_alive_global_start']
        number_alive_global_end = s['number_alive_global_end']
        ENV_START = s['ENV_START']
        RANGE_R = s['RANGE_R']
        ALIVE_THREASHOLD = s['ALIVE_THRESHOLD']
        number_alive_diff.append((int(ENV_START[0]),
                                  int((number_alive_global_start - number_alive_global_end)),
                                  ALIVE_THREASHOLD
                                  ))
        if(ALIVE_THREASHOLD not in alive_threshholds):
            alive_threshholds.append(ALIVE_THREASHOLD)

    finally:
        s.close()

if __name__ == '__main__':

    fig, ax = plt.subplots(figsize=(12, 6))

    plt.title('Species Alive Difference at : ' + str(ALIVE_THREASHOLD))

    plt.xlabel('Environment Start')
    plt.ylabel('Alive Species Difference Start End')
    plt.xlim([-1,RANGE_R+1])


    # - have to lock down to one set of mu and u for alive_threshold diffs

    print(alive_threshholds)

    #print(number_alive_diff)
    #number_alive_diff.sort()
    #print(number_alive_diff)

    for each_threshold in alive_threshholds:
        exp_threshold = []
        x=[]
        y=[]
        for item in number_alive_diff:
            if(item[2] == each_threshold):
                exp_threshold.append(item)

        print(exp_threshold)
        exp_threshold.sort()
        print(exp_threshold)

        for each_entry in exp_threshold:
            x.append(each_entry[0])
            y.append(each_entry[1])

        print(x)
        print(y)
        print("===")
        ax.plot(x,y, label = str(each_threshold))

    plt.legend()
    plt.savefig("3d_alives_diff_" + str(random.randint(100, 999)) + ".png")
    plt.show()

    print("Completed")
