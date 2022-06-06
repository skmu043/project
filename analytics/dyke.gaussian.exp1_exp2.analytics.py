import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

count = 0

for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.gaussian.exp1_exp2.data")
    #print(count)
    try:
        omega = s['omega']
        mu = s['mu']
        niche = s['NICHE']
        survival_threshold = s['SURVIVAL_THRESHOLD']
        env_start = s['ENV_START']
        env_end = s['ENV_END']
        num_alive_start = ['NUMBER_ALIVE_START']
        num_alive_end = s['NUMBER_ALIVE_END']

        RESULT_DATA.append((omega,mu,niche,survival_threshold,env_start,env_end,num_alive_start,num_alive_end))
        count += 1

    finally:
        s.close()

for data_point in RESULT_DATA:
    if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
        UNIQ_SAMPLES.append((data_point[0],data_point[1]))

for sample in UNIQ_SAMPLES:
    print(sample)
#data verification

DATA_VERIFICATION = []

for uniq_ in UNIQ_SAMPLES:
    start_temp_stats = []
    niche_stats = []
    survival_threshold_stats = []

    for data_point in RESULT_DATA:
        if uniq_ == ((data_point[0], data_point[1])):
            if([data_point[2],0] not in niche_stats):
                niche_stats.append([data_point[2], 0])
            if([data_point[3],0] not in survival_threshold_stats):
                survival_threshold_stats.append([data_point[3], 0])
            if([data_point[4],0] not in start_temp_stats):
                start_temp_stats.append([data_point[4], 0])

    start_temp_stats.sort()
    niche_stats.sort()
    survival_threshold_stats.sort()

    index_ = 0
    for start_t in (start_temp_stats):
        for data_point in RESULT_DATA:
            if (uniq_ == ((data_point[0], data_point[1])) and data_point[4] == start_t[0]):
                start_temp_stats[index_][1] += 1

        index_ +=1

    index_ = 0
    for niche_s in (niche_stats):
        for data_point in RESULT_DATA:
            if (uniq_ == ((data_point[0], data_point[1])) and data_point[2] == niche_s[0]):
                niche_stats[index_][1] += 1

        index_ +=1

    index_ = 0
    for st_v in (survival_threshold_stats):
        for data_point in RESULT_DATA:
            if (uniq_ == ((data_point[0], data_point[1])) and data_point[3] == st_v[0]):
                survival_threshold_stats[index_][1] += 1

        index_ +=1

    #print(uniq_)
    #print(start_temp_stats)
    #print(niche_stats)
    #print(survival_threshold_stats)

    DATA_VERIFICATION.append([start_temp_stats,niche_stats,survival_threshold_stats])

temp_check = []
niche_check = []
st_check = []

#for sample in DATA_VERIFICATION:
#    for each_temp in sample[0]:
#        temp_check.append(each_temp[1])
#    for each_niche in sample[1]:
#        niche_check.append(each_niche[1])
#    for each_st in sample[2]:
#        st_check.append(each_st[1])

#print(all(i == temp_check[0] for i in temp_check))
#print(all(i == niche_check[0] for i in niche_check))
#print(all(i == st_check[0] for i in st_check))

print("DATA VALIDATION CHECK : ")
print(all(i == DATA_VERIFICATION[0] for i in DATA_VERIFICATION))

print("=====================")

XFONT = 40
YFONT = 40
X_TICKS = 20
Y_TICKS = 20
XFIG = 30
YFIG = 30
TFONT = 40


def end_temp():
    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('START/END Temp : 0.2 start temp 50', fontsize=TFONT)
    plt.xlabel('With Truncation', fontsize=XFONT)
    plt.ylabel('Zero Truncation', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.ylim(-50, 150)
    plt.xlim(-50, 150)

    for data_point in RESULT_DATA:
        x = 0
        y = 0

        if(data_point[3] == 0 and data_point[4] == 50 and data_point[2] == 5):
            x = data_point[5]
            for data_point_2 in RESULT_DATA:
                if(data_point[0] == data_point_2[0]
                        and data_point[1] == data_point_2[1]
                        and data_point_2[3] == 0.2
                        and data_point_2[4] == 50
                        and data_point_2[2] == 5):
                    y = data_point_2[5]

        plt.scatter(x,y)

    plt.show()


#=======================================================================================================================
def number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end():

    x = [] #start temp
    y = [] #number alive at the end
    z = [] #survival threshold

    #RESULT_DATA.append((omega,mu,niche,survival_threshold,env_start,env_end,num_alive_start,num_alive_end))

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            x.append(data_point[4])
            y.append(data_point[7])
            z.append(data_point[3])

    uniq_start_temps = np.unique(np.array(x))
    uniq_survivals = np.unique(np.array(z))
    main_result = []

    print(uniq_start_temps)
    print(uniq_survivals)

    for each_survival_threshold in uniq_survivals:

        plt.figure(figsize=(XFIG,YFIG), dpi=200)
        plt.title('The number of alive species at the end of simulations: ' + str(each_survival_threshold), fontsize=TFONT)
        plt.xlabel('Starting Temperature', fontsize=XFONT)
        plt.ylabel('Number Alive at the End', fontsize=YFONT)
        plt.xticks(fontsize=X_TICKS)
        plt.yticks(fontsize=Y_TICKS)

        alive_below = []
        alive_above = []

        for each_start_temp in uniq_start_temps:
            index_1 = 0
            below_zero = 0
            above_zero = 0
            for each_row in x:
                if(x[index_1] == each_start_temp and z[index_1] == each_survival_threshold):
                    if(y[index_1] <= 0):
                        below_zero +=1
                    if(y[index_1] > 0):
                        above_zero +=1
                index_1 +=1

            alive_below.append(below_zero)
            alive_above.append(above_zero)

            main_result.append([each_start_temp,each_survival_threshold,below_zero,above_zero])

        for each in main_result:
            print(each)

        X=[]
        for each in uniq_start_temps:
            X.append(str(each))

        X_axis = np.arange(len(X))

        plt.bar(X_axis - 0.2, alive_below, 0.4, label = 'Simulations with no alives at the end')
        plt.bar(X_axis + 0.2, alive_above, 0.4, label = 'Simulation with alives at the end')
        plt.xticks(X_axis, X)

        plt.legend()
        plt.show()
        #[uniq_start_temp (19), survival_threshold (2),number of simulations with zero alives at end, number of alives above 0]


#=======================================================================================================================
number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end()
#=======================================================================================================================
def zero_alives_at_end_removed_number_alive_at_each_start_temperature_at_the_end_of_simulation():
    print()
