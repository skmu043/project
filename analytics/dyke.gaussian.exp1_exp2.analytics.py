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
        num_alive_start = s['NUMBER_ALIVE_START']
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
#number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end()
#=======================================================================================================================
def number_alive_at_each_start_temperature_at_the_start_of_simulation():


    # Write UP notes - if simulation ended with zero species alive - its still included in these results
    # Where it won't be included is in the temperature section where species alive at the end matters
    # no species alive at the end of a simulation and a temperature between 0 and R does not mean the
    # simulation is good - there cannot be regulation without species (temperature stays at a value with no change
    # when the species no longer exist)

    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))

    start_temp_0 = []
    alive_start_0 = []
    start_temp_2 = []
    alive_start_2 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp_0.append(data_point[4])
                alive_start_0.append(data_point[6])
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                alive_start_2.append(data_point[6])


    print(start_temp_0)
    print(alive_start_0)
    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('The Number of Alive at the Start of the Simulation',fontsize=30)
    #fig.set_size_inches(3, 1.5)
    ax1.scatter(start_temp_0, alive_start_0)
    ax1.set_title('The Dyke Weaver Model', fontsize=20)
    ax1.set_xlabel('Temperature', fontsize=20)
    ax1.set_ylabel('Number Alives', fontsize=20)
    ax2.scatter(start_temp_2, alive_start_2)
    ax2.set_title('Survival Threshold of 0.2', fontsize=20)
    ax2.set_xlabel('Temperature', fontsize=20)
    ax2.set_ylabel('Number Alives', fontsize=20)
    ax2.legend(prop={'size': 20})
    fig.show()
#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
def number_alive_at_each_start_temperature_at_the_end_of_simulation():


    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))

    start_temp_0 = []
    alive_end_0 = []
    start_temp_2 = []
    alive_end_2 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp_0.append(data_point[4])
                alive_end_0.append(data_point[7])
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                alive_end_2.append(data_point[7])


    print(start_temp_0)
    print(alive_end_0)
    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('The Number of Alive at the End of the Simulation',fontsize=30)
    #fig.set_size_inches(3, 1.5)
    ax1.scatter(start_temp_0, alive_end_0)
    ax1.set_title('The Dyke Weaver Model', fontsize=20)
    ax1.set_xlabel('Temperature', fontsize=20)
    ax1.set_ylabel('Number Alives', fontsize=20)
    ax2.scatter(start_temp_2, alive_end_2)
    ax2.set_title('Survival Threshold of 0.2', fontsize=20)
    ax2.set_xlabel('Temperature', fontsize=20)
    ax2.set_ylabel('Number Alives', fontsize=20)
    ax2.legend(prop={'size': 20})
    fig.show()
#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================
def average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation():

    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))

    start_temp_0 = []
    alives_start_0 = []
    alive_end_0 = []
    start_temp_2 = []
    alives_start_2 = []
    alive_end_2 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp_0.append(data_point[4])
                alives_start_0.append(data_point[6])
                alive_end_0.append(data_point[7])
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                alives_start_2.append(data_point[6])
                alive_end_2.append(data_point[7])


    print(start_temp_0)
    print(alive_end_0)
    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('Average Number of alives at the Start and End of the Simulation',fontsize=30)
    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(start_temp_0))
    uniq_start_temps.sort()

    print(uniq_start_temps)


    # alive end stats
    stats_temp = []
    stats_avg_alives = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += alive_end_0[index]
                count += 1

            index += 1
        stats_temp.append(each_temp)
        stats_avg_alives.append((sum/count))

    ax1.plot(stats_temp, stats_avg_alives, label='alives end')

    stats_temp = []
    stats_avg_alives = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += alives_start_0[index]
                count += 1

            index += 1
        stats_temp.append(each_temp)
        stats_avg_alives.append((sum/count))

    ax1.plot(stats_temp, stats_avg_alives, label='alives start')



    ax1.set_title('The Dyke Weaver Model', fontsize=20)
    ax1.set_xlabel('Temperature', fontsize=20)
    ax1.set_ylabel('Number Alives', fontsize=20)

    stats_temp_2 = []
    stats_avg_alives_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += alive_end_2[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_alives_2.append((sum/count))


    ax2.plot(stats_temp_2, stats_avg_alives_2,label='alives end')

    stats_temp_2 = []
    stats_avg_alives_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += alives_start_2[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_alives_2.append((sum/count))


    ax2.plot(stats_temp_2, stats_avg_alives_2, label='alives start')
    ax2.set_title('Survival Threshold of 0.2', fontsize=20)
    ax2.set_xlabel('Temperature', fontsize=20)
    ax2.set_ylabel('Number Alives', fontsize=20)
    ax2.legend(prop={'size': 20})
    fig.show()
#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation()
#=======================================================================================================================

def average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels():

    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))

    start_temp_0 = []
    alives_start_0 = []
    alive_end_0 = []
    survival_thresh_0 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5):
            start_temp_0.append(data_point[4])
            survival_thresh_0.append(data_point[3])
            alives_start_0.append(data_point[6])
            alive_end_0.append(data_point[7])

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, dpi=300, figsize=(30,10))
    fig.suptitle('Average Number of alives at the Start and End of the Simulation for Survival Thresholds',fontsize=30)
    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(start_temp_0))
    uniq_start_temps.sort()
    survival_thresh_0.sort()

    print(uniq_start_temps)
    print(survival_thresh_0)

    # alive end stats

    st2_s = []
    st4_s = []
    st6_s = []
    st8_s = []

    st2_e = []
    st4_e = []
    st6_e = []
    st8_e = []

    for each_st in survival_thresh_0:
        stats_temp = []
        stats_avg_alives = []
        for each_temp in uniq_start_temps:
            index = 0
            sum = 0
            count = 0
            for each_instance in start_temp_0:
                if(each_temp == each_instance and each_st == survival_thresh_0[index]):
                    sum += alive_end_0[index]
                    count += 1
                index += 1
            stats_temp.append(each_temp)
            stats_avg_alives.append((sum/count))

        if(each_st == 0.2):
            st2_e = stats_avg_alives.copy()
            ax1.plot(stats_temp, stats_avg_alives, label='alives end 0.2')
        if(each_st == 0.4):
            st4_e = stats_avg_alives.copy()
            ax2.plot(stats_temp, stats_avg_alives, label='alives end 0.4')
        if(each_st == 0.6):
            st6_e = stats_avg_alives.copy()
            ax3.plot(stats_temp, stats_avg_alives, label='alives end 0.6')
        if(each_st == 0.8):
            st8_e = stats_avg_alives.copy()
            ax4.plot(stats_temp, stats_avg_alives, label='alives end 0.8')

        stats_temp = []
        stats_avg_alives = []

        for each_temp in uniq_start_temps:
            index = 0
            sum = 0
            count = 0
            for each_instance in start_temp_0:
                if(each_temp == each_instance and each_st == survival_thresh_0[index]):
                    sum += alives_start_0[index]
                    count += 1
                index += 1
            stats_temp.append(each_temp)
            stats_avg_alives.append((sum/count))

        if(each_st == 0.2):
            st2_s = stats_avg_alives.copy()
            ax1.plot(stats_temp, stats_avg_alives, label='alives start 0.2')
            ax1.set_title('T0.2', fontsize=20)
            ax1.set_xlabel('Temperature', fontsize=20)
            ax1.set_ylabel('Number Alives', fontsize=20)
        if(each_st == 0.4):
            st4_s = stats_avg_alives.copy()
            ax2.plot(stats_temp, stats_avg_alives, label='alives start 0.4')
            ax2.set_title('T0.4', fontsize=20)
            ax2.set_xlabel('Temperature', fontsize=20)
            ax2.set_ylabel('Number Alives', fontsize=20)
        if(each_st == 0.6):
            st6_s = stats_avg_alives.copy()
            ax3.plot(stats_temp, stats_avg_alives, label='alives start 0.6')
            ax3.set_title('T0.6', fontsize=20)
            ax3.set_xlabel('Temperature', fontsize=20)
            ax3.set_ylabel('Number Alives', fontsize=20)
        if(each_st == 0.8):
            st8_s = stats_avg_alives.copy()
            ax4.plot(stats_temp, stats_avg_alives, label='alives start 0.8')
            ax4.set_title('T0.8', fontsize=20)
            ax4.set_xlabel('Temperature', fontsize=20)
            ax4.set_ylabel('Number Alives', fontsize=20)

    fig.show()


    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('Number Alive for each Temperature at each Survival Threshold',fontsize=30)
    #fig.set_size_inches(3, 1.5)

    index = 0
    for temp in uniq_start_temps:
        linex = [0.2,0.4,0.6,0.8]
        liney = [st2_s[index],st4_s[index],st6_s[index],st8_s[index]]
        ax2.plot(linex, liney, label=str(temp))
        index+=1


    ax1.set_title('Alives Start', fontsize=20)
    ax1.set_xlabel('Survival Threshold', fontsize=20)
    ax1.set_ylabel('Number Alives', fontsize=20)

    index = 0
    for temp in uniq_start_temps:
        linex = [0.2,0.4,0.6,0.8]
        liney = [st2_e[index],st4_e[index],st6_e[index],st8_e[index]]
        ax1.plot(linex, liney, label=str(temp))
        index+=1

    ax2.set_title('Alives End', fontsize=20)
    ax2.set_xlabel('Survival Threshold', fontsize=20)
    ax2.set_ylabel('Number Alives', fontsize=20)
    ax2.legend(prop={'size': 20})
    fig.show()

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Different Start and End Alives at Suvival Thresholds', fontsize=TFONT)
    plt.xlabel('Survival Threshold', fontsize=XFONT)
    plt.ylabel('Number Alive', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    index = 0
    for temp in uniq_start_temps:
        linex = [0.2,0.4,0.6,0.8]
        liney = [st2_e[index]-st2_s[index],st4_e[index]-st4_s[index],st6_e[index]-st6_s[index],st8_e[index]-st8_s[index]]
        plt.plot(linex, liney, label=str(temp))
        index+=1
    plt.legend()
    plt.show()


#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels()
#=======================================================================================================================

#=======================================================================================================================
def number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R():

    x = [] #start temp
    y = [] #env_end
    z = [] #survival threshold
    al = []

    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            x.append(data_point[4])
            y.append(data_point[5])
            z.append(data_point[3])
            al.append(data_point[7])

    uniq_start_temps = np.unique(np.array(x))
    uniq_start_temps.sort()
    uniq_survivals = np.unique(np.array(z))
    uniq_survivals.sort()

    main_result = []

    print(uniq_start_temps)
    print(uniq_survivals)

    for each_survival_threshold in uniq_survivals:

        plt.figure(figsize=(XFIG,YFIG), dpi=200)
        plt.title('Bounds [inside 0-R or outside 0-R] : ' + str(each_survival_threshold), fontsize=TFONT)
        plt.xlabel('Start Temperature', fontsize=XFONT)
        plt.ylabel('Simulations Bounds', fontsize=YFONT)
        plt.xticks(fontsize=X_TICKS)
        plt.yticks(fontsize=Y_TICKS)

        inside_bounds = []
        outside_bounds = []

        for each_start_temp in uniq_start_temps:
            index_1 = 0
            inside_b = 0
            outside_b = 0
            for each_row in x:
                if(x[index_1] == each_start_temp and z[index_1] == each_survival_threshold and al[index_1] > 0): # al = Number of alive species greater than one
                    if(y[index_1] > 0 and y[index_1] < 100):
                        inside_b +=1
                    if(y[index_1] < 0 or y[index_1] > 100):
                        outside_b +=1
                index_1 +=1

            inside_bounds.append(inside_b)
            outside_bounds.append(outside_b)

            main_result.append([each_start_temp,each_survival_threshold,inside_b,outside_b])

        X=[]
        for each in uniq_start_temps:
            X.append(str(each))

        X_axis = np.arange(len(X))

        plt.bar(uniq_start_temps - 0.5, inside_bounds, 1, label = 'Simulations within bounds 0 - R')
        plt.bar(uniq_start_temps + 0.5, outside_bounds, 1, label = 'Simulations outside bounds 0 - R')


        #plt.xticks(X_axis, X)

        plt.legend(prop={'size': 20})
        plt.show()
        #[uniq_start_temp (19), survival_threshold (2),number of simulations with zero alives at end, number of alives above 0]


#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R()
#=======================================================================================================================

#=======================================================================================================================
def number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only():

    UNIQ_SAMPLES = []
    for data_point in RESULT_DATA:
        if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
            UNIQ_SAMPLES.append((data_point[0],data_point[1]))

    mu = []
    omega = []
    x = [] #start temp
    y = [] #env_end
    z = [] #survival threshold
    al = []

    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            omega.append(data_point[0])
            mu.append(data_point[1])
            x.append(data_point[4])  #env_start
            y.append(data_point[5])  #env_end
            z.append(data_point[3])  #survival_threshold
            al.append(data_point[7]) #num_alive_end

    uniq_start_temps = np.unique(np.array(x))
    uniq_start_temps.sort()
    uniq_survivals = np.unique(np.array(z))
    uniq_survivals.sort()

    main_result = []

    print(uniq_start_temps)
    print(uniq_survivals)

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Bounds [inside 0-R or outside 0-R] :', fontsize=TFONT)
    plt.xlabel('Start Temperature', fontsize=XFONT)
    plt.ylabel('Simulations Bounds', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    main_results = []

    for each_temp in uniq_start_temps:
        both_all = []
        zero_all = []
        zero2_all = []
        for each_sample in UNIQ_SAMPLES:
            index = 0
            both = 0
            zero = 0
            zero2 = 0
            for each_data in x:
                if(x[index] == each_temp and ((omega[index], mu[index])) == each_sample and al[index] > 0):
                    if(z[index]==0 and (y[index] > 0 and y[index] < 100)):
                        zero +=1
                    if(z[index]==0.2 and (y[index] > 0 and y[index] < 100)):
                        zero2 +=1
                index +=1
            if(zero > 0 and zero2 > 0):
                both_all.append(1)
            if(zero > 0 and zero2 == 0):
                zero_all.append(1)
            if(zero == 0 and zero2 > 0):
                zero2_all.append(1)

        main_results.append((each_temp,sum(both_all), sum(zero_all), sum(zero2_all)))

    bar_both = []
    bar_zero_all = []
    bar_zero2_all = []

    for each_result in main_results:
        bar_both.append(each_result[1])
        bar_zero_all.append(each_result[2])
        bar_zero2_all.append(each_result[3])


    plt.bar(uniq_start_temps - 0.7, bar_both, 1, label = 'Both within bounds 0 - R')
    plt.bar(uniq_start_temps + 0.5, bar_zero_all, 1, label = 'zero only inside bounds 0 - R')
    plt.bar(uniq_start_temps + 0.7, bar_zero2_all, 1, label = '0.2 only inside bounds 0 - R')
    plt.legend(prop={'size': 20})
    plt.show()


#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only()
#=======================================================================================================================

def end_temperature_both_zero_zero2_with_alives_overlay():

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('START/END Temp : 0.2 start temp 50', fontsize=TFONT)
    plt.xlabel('With Truncation', fontsize=XFONT)
    plt.ylabel('Zero Truncation', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.ylim(-50, 150)
    plt.xlim(-50, 150)


    #RESULT_DATA.append((
    # omega[0],
    # mu[1],
    # niche[2],
    # survival_threshold[3],
    # env_start[4],
    # env_end[5],
    # num_alive_start[6],
    # num_alive_end[7]
    # ))
    zero_t = []
    two_t = []
    alives = []

    for data_point in RESULT_DATA:

        if((data_point[3] == 0 or data_point[3] == 0.2) and data_point[4] == 50 and data_point[2] == 5):
            if(data_point[3] == 0):
                zero_t.append(data_point[5])
            if(data_point[3] == 0.2):
                two_t.append(data_point[5])
                alives.append(data_point[7])

    plt.scatter(zero_t, two_t, c=alives, cmap = 'viridis')
    plt.colorbar(label="Number Alive", orientation="vertical")
    plt.show()
#=======================================================================================================================
end_temperature_both_zero_zero2_with_alives_overlay()
#=======================================================================================================================

