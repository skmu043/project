import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

number_alive_global_start_ = []
number_alive_local_start_ = []
number_alive_global_end_ = []
number_alive_local_end_ = []
number_alive_start_ = []
number_alive_end_ = []
ENV_START_ = []

for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor.rk4.data")

    try:

        #RUN_ID = s['RUN_ID']

        number_alive_global_start = s['number_alive_global_start']
        number_alive_local_start = s['number_alive_local_start']
        number_alive_global_end = s['number_alive_global_end']
        number_alive_local_end = s['number_alive_local_end']

        number_alive_start = s['number_alive_start']
        number_alive_end = s['number_alive_end']

        ENV_START = s['ENV_START']

        RANGE_R = s['RANGE_R']


        number_alive_global_start_.append(number_alive_global_start)
        number_alive_local_start_.append(number_alive_local_start)
        number_alive_global_end_.append(number_alive_global_end)
        number_alive_local_end_.append(number_alive_local_end)
        number_alive_start_.append(number_alive_start)
        number_alive_end_.append(number_alive_end)
        ENV_START_.append(ENV_START)

        # SUPER DATA STRUCTURE NEEDED

        #(Eg, El, number_alive_global_start, number_alive_local_start, number_alive_start, number_alive_global_end, number_alive_local_end, number_alive_end)

        #Eg = [] ...
        #El = [] ...
        #number_alive_global_start= []

        # Loops would be for item in Eg ...
        # plot Eg = x , El = y (the correct one) , heatmap += number_alive_global_start

        # Replicate for the ones below :

        #number_alive_local_start
        #number_alive_start
        #number_alive_global_end
        #number_alive_local_end
        #number_alive_end

    finally:
        s.close()

if __name__ == '__main__':

    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Heat Map of Alives for Global Start', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)

    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]

    idx = 0
    for start_vars in ENV_START_:
        heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_global_start_[idx]
        idx += 1


    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')

    plt.show()


    print("Completed")
