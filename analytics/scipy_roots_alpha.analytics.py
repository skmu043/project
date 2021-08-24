import shelve, os
import matplotlib.pyplot as plt
import sys
import numpy as np
print(sys.version)
import math, random, time
import statistics


K       = 0
N       = 0
RUN_ID  = 0

data_points = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    s = shelve.open(data_dr + "/" + str(si) + "/scipy_roots_alpha.data")
    try :
        args            = s['sys.argv']
        w               = s['w']
        u               = s['u']
        K               = s['K']
        N               = s['N']
        stable_point    = s['stable_point']
        RUN_ID          = s['RUN_ID']

        data_points.append((K, N, stable_point))

    finally:
        s.close()

uniq_n = []
for x in data_points:
    if x[1] not in uniq_n:
        uniq_n.append(x[1])

uniq_n.sort()
print(uniq_n)
uniq_k = []
for x in data_points:
    if x[0] not in uniq_k:
        uniq_k.append(x[0])

uniq_k.sort()
print(uniq_k)

def plot_stable_biotic():
    ax = plt.figure(figsize=(20,10))
    plt.title('Number of Species vs Stable Points', fontsize=40)
    plt.xlabel('Species', fontsize=20)
    plt.ylabel('Stable Points (mean)', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.axvline(x=0)
    plt.axhline(y=0)
    plt.xticks([2, 5, 10, 20, 50, 100, 200])

    # This Set Does Total Abundance Sum, Mean and Standard Deviation (Error Bar Plot)
    x = np.array([])
    y = np.array([])
    e = np.array([])

    biotic_elements_x   = []
    stable_points_y     = []

    #for each_N in uniq_n:
    #    for each_data_point in data_points:
    #        if(each_data_point[1] == each_N):
    #            biotic_elements_x.append(each_data_point[0])
    #            stable_points_y.append(each_data_point[2])
        #plt.plot(biotic_elements_x,stable_points_y, '.',label = each_N)

    #    biotic_elements_x.clear()
    #    stable_points_y.clear()

    #print(data_points)
    stable_points_data = []

    k_x = []
    s_y = []
    plot_k_x = []
    plot_s_y = []

    for each_N in uniq_n:
        for each_K in uniq_k:
            for each_data_point in data_points:
                # (K, N, stable_point)
                if(each_K == each_data_point[0] and each_N == each_data_point[1]):
                    stable_points_data.append(each_data_point[2])
            k_x.append(each_K)
            s_y.append(statistics.mean(stable_points_data))
            stable_points_data.clear()

            #x = np.append(x, each_K)
            #y = np.append(y, statistics.mean(stable_points_data))
            #e = np.append(e, statistics.stdev(stable_points_data))
            #print("N: ", each_N)
            #print("K: ", k_x)
            #print("S: ", s_y)
            #print()
        plot_k_x.append(k_x.copy())
        plot_s_y.append(s_y.copy())
        k_x.clear()
        s_y.clear()
        index_num = 0
        #print("PlotK : ", plot_k_x)
        #print("PlotS : ", plot_s_y)
    marker = ""
    for each_N in uniq_n:
        if(each_N == 5):
            marker = "r-"
        elif(each_N == 7):
            marker = "b-"
        elif(each_N == 10):
            marker = "k-"
        plt.plot(plot_k_x[index_num], plot_s_y[index_num], marker, label = str(each_N))
        index_num +=1
        #plt.errorbar(x, y, e, linestyle='None', marker='^', elinewidth=7, capsize=8, capthick=7)


    plt.legend(loc="upper right")
    plt.savefig('stable_biotic_x.png')
    plt.show()

plot_stable_biotic()


