import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm

RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

count = 0

for file in tqdm(data_archives):
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
        total_abundance_start = s['TOTAL_ABUNDANCE_START']
        total_abundance_end = s['TOTAL_ABUNDANCE_END']

        RESULT_DATA.append((omega,mu,niche,survival_threshold,env_start,env_end,num_alive_start,num_alive_end,total_abundance_start,total_abundance_end))
        count += 1

    finally:
        s.close()

#=======================================================================================================================
def data_verification():

    for data_point in RESULT_DATA:
        if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
            UNIQ_SAMPLES.append((data_point[0],data_point[1]))

    #for sample in UNIQ_SAMPLES:
    #    print(sample)
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
#data_verification()
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

    #print(uniq_start_temps)
    #print(uniq_survivals)

    for each_survival_threshold in uniq_survivals:

        plt.figure(figsize=(XFIG,YFIG), dpi=200)
        plt.title('Dyke/Weaver Model : Simulations with Alive Species', fontsize=TFONT)
        if(each_survival_threshold == 0.2):
            plt.title('Simulations with Alive Species using a survival threshold of 0.2', fontsize=TFONT)
        plt.xlabel('Starting Temperature', fontsize=XFONT)
        plt.ylabel('Simulations with/without Alive Species', fontsize=YFONT)
        plt.xticks(fontsize=X_TICKS)
        plt.yticks(fontsize=Y_TICKS)

        alive_below = []
        alive_above = []

        for each_start_temp in uniq_start_temps:
            index_1 = 0
            below_zero = 0
            above_zero = 0
            for each_row in x:
                if(each_row == each_start_temp and z[index_1] == each_survival_threshold):
                    if(y[index_1] <= 0):
                        below_zero +=1
                    if(y[index_1] > 0):
                        above_zero +=1
                    #print(each_survival_threshold,each_start_temp,x[index_1],y[index_1],z[index_1])
                index_1 +=1

            alive_below.append(below_zero)
            alive_above.append(above_zero)

            #main_result.append([each_start_temp,each_survival_threshold,below_zero,above_zero])

        X=[]
        for each in uniq_start_temps:
            X.append(str(each))

        X_axis = np.arange(len(X))

        #print(main_result)

        plt.bar(X_axis - 0.2, alive_below, 0.4, label = 'Simulations with no alive species')
        plt.bar(X_axis + 0.2, alive_above, 0.4, label = 'Simulation with alive species')
        plt.xticks(X_axis, X)

        plt.legend(prop={'size': 30})
        plt.tight_layout()
        plt.savefig('number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end_'+str(each_survival_threshold)+'.jpg')
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


    #print(start_temp_0)
    #print(alive_start_0)

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Number of alive species at the start of the simulation', fontsize=TFONT)
    plt.xlabel('Starting Temperature', fontsize=XFONT)
    plt.ylabel('Alive Species', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.scatter(start_temp_0, alive_start_0, label='Dyke/Weaver Model')
    plt.scatter(start_temp_2, alive_start_2, label="Survival Threshold 0.2")
    plt.xticks(np.arange(0,100,5))
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_start_of_simulation.jpg')
    plt.show()
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


    #print(start_temp_0)
    #print(alive_end_0)

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Number of alive species at the end of the simulation', fontsize=TFONT)
    plt.xlabel('Starting Temperature', fontsize=XFONT)
    plt.ylabel('Alive Species', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.scatter(start_temp_0, alive_end_0, label='Dyke/Weaver Model')
    plt.scatter(start_temp_2, alive_end_2, label="Survival Threshold 0.2")
    plt.xticks(np.arange(0,100,5))
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_end_of_simulation.jpg')
    plt.show()


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


    #print(start_temp_0)
    #print(alive_end_0)

    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(start_temp_0))
    uniq_start_temps.sort()

    #print(uniq_start_temps)


    # alive end stats
    stats_temp = []
    stats_avg_alives = []

    dw_temp_s = []
    dw_alive_s = []
    dw_temp_e = []
    dw_alive_e = []

    st_temp_s = []
    st_alive_s = []
    st_temp_e = []
    st_alive_e = []


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

    dw_temp_e = stats_temp.copy()
    dw_alive_e = stats_avg_alives.copy()

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

    dw_temp_s = stats_temp.copy()
    dw_alive_s = stats_avg_alives.copy()

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

    st_temp_e = stats_temp_2.copy()
    st_alive_e = stats_avg_alives_2.copy()

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

    st_temp_s = stats_temp_2.copy()
    st_alive_s = stats_avg_alives_2.copy()



    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Average Number of alives at the Start and End of the Simulation', fontsize=TFONT)
    plt.xlabel('Starting Temperature', fontsize=XFONT)
    plt.ylabel('Average Number of Alive Species', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.plot(dw_temp_s, dw_alive_s, label='Dyke/Weaver Alive Species at Start')
    plt.plot(dw_temp_e, dw_alive_e, label='Dyke/Weaver Alive Species at End')
    plt.plot(st_temp_s, st_alive_s, label='Survival Threshold 0.2 Alive Species at Start')
    plt.plot(st_temp_e, st_alive_e, label='Survival Threshold 0.2 Alive Species at End')
    plt.xticks(np.arange(0,100,5))
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation.jpg')
    plt.show()


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

    fig = plt.figure(figsize=(XFIG,YFIG),dpi=200)
    spec = fig.add_gridspec(ncols=2, nrows=2)
    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[1, 0])
    ax4 = fig.add_subplot(spec[1, 1])
    fig.suptitle('Average Number of alive species at the Start and End of simulations for different Survival Thresholds',fontsize=30)

    uniq_start_temps = np.unique(np.array(start_temp_0))
    uniq_start_temps.sort()
    survival_thresh_0.sort()

    #print(uniq_start_temps)
    #print(survival_thresh_0)

    # alive end stats

    st2_s = []
    st4_s = []
    st6_s = []
    st8_s = []

    st2_e = []
    st4_e = []
    st6_e = []
    st8_e = []

    #print("start")
    for each_st in tqdm(survival_thresh_0):
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
            ax1.plot(stats_temp, stats_avg_alives, label='Alive Species at End')
        if(each_st == 0.4):
            st4_e = stats_avg_alives.copy()
            ax2.plot(stats_temp, stats_avg_alives, label='Alive Species at End')
        if(each_st == 0.6):
            st6_e = stats_avg_alives.copy()
            ax3.plot(stats_temp, stats_avg_alives, label='Alive Species at End')
        if(each_st == 0.8):
            st8_e = stats_avg_alives.copy()
            ax4.plot(stats_temp, stats_avg_alives, label='Alive Species at End')

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
            ax1.plot(stats_temp, stats_avg_alives, label='Alive Species at Start')
            ax1.set_title('Survival Threshold 0.2',fontsize=25)
            ax1.set_xlabel('Temperature',fontsize=15)
            ax1.set_ylabel('Average number of alive species',fontsize=15)
        if(each_st == 0.4):
            st4_s = stats_avg_alives.copy()
            ax2.plot(stats_temp, stats_avg_alives, label='Alive Species at Start')
            ax2.set_title('Survival Threshold 0.4',fontsize=25)
            ax2.set_xlabel('Temperature',fontsize=15)
            ax2.set_ylabel('Average number of alive species',fontsize=15)
        if(each_st == 0.6):
            st6_s = stats_avg_alives.copy()
            ax3.plot(stats_temp, stats_avg_alives, label='Alive Species at Start')
            ax3.set_title('Survival Threshold 0.6',fontsize=25)
            ax3.set_xlabel('Temperature',fontsize=15)
            ax3.set_ylabel('Average number of alive species',fontsize=15)
        if(each_st == 0.8):
            st8_s = stats_avg_alives.copy()
            ax4.plot(stats_temp, stats_avg_alives, label='Alive Species at Start')
            ax4.set_title('Survival Threshold 0.8',fontsize=25)
            ax4.set_xlabel('Temperature',fontsize=15)
            ax4.set_ylabel('Average number of alive species',fontsize=15)

    fig.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels_quad.jpg')
    fig.show()

    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('Number Alive for each Temperature at each Survival Threshold',fontsize=30)
    #fig.set_size_inches(3, 1.5)

    #plt.bar(X_axis - 0.2, alive_below, 0.4, label = 'Simulations with no alive species')


    index = 0
    for temp in uniq_start_temps:
        linex = [0.2,0.4,0.6,0.8]
        liney = [st2_s[index],st4_s[index],st6_s[index],st8_s[index]]
        ax2.scatter(linex, liney, label=str(temp))
        index+=1


    ax1.set_title('Species at the start of simulations', fontsize=20)
    ax1.set_xlabel('Survival Threshold', fontsize=20)
    ax1.set_ylabel('Alive Species', fontsize=20)

    index = 0
    for temp in uniq_start_temps:
        linex = [0.2,0.4,0.6,0.8]
        liney = [st2_e[index],st4_e[index],st6_e[index],st8_e[index]]
        ax1.scatter(linex, liney, label=str(temp))
        index+=1

    ax2.set_title('Species at the end of simulations', fontsize=20)
    ax2.set_xlabel('Survival Threshold', fontsize=20)
    ax2.set_ylabel('Alive Species', fontsize=20)
    ax2.legend()
    fig.tight_layout()
    fig.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels_threshold.jpg')
    fig.show()

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Difference between the start and end number of alive species', fontsize=TFONT)
    plt.xlabel('Survival Threshold', fontsize=XFONT)
    plt.ylabel('Alive Species', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    index = 0
    for temp in uniq_start_temps:
        linex = [0.2,0.4,0.6,0.8]
        liney = [st2_e[index]-st2_s[index],st4_e[index]-st4_s[index],st6_e[index]-st6_s[index],st8_e[index]-st8_s[index]]
        plt.scatter(linex, liney, label=str(temp))
        index+=1
    plt.legend()
    plt.tight_layout()
    plt.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels.jpg')
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

    for each_survival_threshold in tqdm(uniq_survivals):

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
        plt.tight_layout()
        plt.savefig('number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R_' + str(each_survival_threshold) + '.jpg' )
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

    for each_temp in tqdm(uniq_start_temps):
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

    width = 1

    plt.bar(uniq_start_temps - width/3, bar_both, width, label = 'Both within bounds 0 - R')
    plt.bar(uniq_start_temps + 0.8, bar_zero_all, width, label = 'Dyke/Weaver only inside bounds 0 - R')
    plt.bar(uniq_start_temps + 2, bar_zero2_all, width, label = 'Survival Threshold only inside bounds 0 - R')
    plt.legend(prop={'size': 25})
    plt.tight_layout()
    plt.savefig('number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only.jpg')
    plt.show()


#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only()
#=======================================================================================================================

def end_temperature_both_zero_zero2_with_alives_overlay():

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('End temperature of Dyke/Weaver vs Survival Threshold at 0.2 starting at 50', fontsize=TFONT)
    plt.xlabel('End Temperature using Survival Threshold', fontsize=XFONT)
    plt.ylabel('End Temperature Dyke/Weaver Model', fontsize=YFONT)
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

    for data_point in tqdm(RESULT_DATA):
        if((data_point[3] == 0 or data_point[3] == 0.2) and data_point[4] == 50 and data_point[2] == 5):
            if(data_point[3] == 0):
                zero_t.append(data_point[5])
            if(data_point[3] == 0.2):
                two_t.append(data_point[5])
                alives.append(data_point[7])

    plt.scatter(zero_t, two_t, c=alives, cmap = 'viridis', s=200)
    cb = plt.colorbar(orientation="vertical")
    cb.ax.tick_params(labelsize=20)
    cb.set_label(label="Number of alive species",size=YFONT)
    plt.tight_layout()
    plt.savefig('end_temperature_both_zero_zero2_with_alives_overlay.jpg')
    plt.show()

#=======================================================================================================================
#end_temperature_both_zero_zero2_with_alives_overlay()
#=======================================================================================================================
def abundance_alive_at_each_start_temperature_at_the_start_of_simulation():


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
    # total_abundance_start[8]
    # total_abundance_end[9])
    # ))

    start_temp_0 = []
    abundance_start_0 = []
    start_temp_2 = []
    abundance_start_2 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp_0.append(data_point[4])
                abundance_start_0.append(data_point[8])
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                abundance_start_2.append(data_point[8])


    #print(start_temp_0)
    #print(alive_start_0)

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Total Abundance of species at the start of the simulation', fontsize=TFONT)
    plt.xlabel('Starting Temperature', fontsize=XFONT)
    plt.ylabel('Total Abundance', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.scatter(start_temp_0, abundance_start_0, label='Dyke/Weaver Model')
    plt.scatter(start_temp_2, abundance_start_2, label="Survival Threshold 0.2")
    plt.xticks(np.arange(0,100,5))
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('abundance_alive_at_each_start_temperature_at_the_start_of_simulation.jpg')
    plt.show()
#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
def abundance_alive_at_each_start_temperature_at_the_end_of_simulation():


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
    # total_abundance_start[8]
    # total_abundance_end[9])
    # ))

    start_temp_0 = []
    abundance_end_0 = []
    start_temp_2 = []
    abundance_end_2 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp_0.append(data_point[4])
                abundance_end_0.append(data_point[9])
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                abundance_end_2.append(data_point[9])


    #print(start_temp_0)
    #print(alive_start_0)

    plt.figure(figsize=(XFIG,YFIG), dpi=200)
    plt.title('Total Abundance of species at the start of the simulation', fontsize=TFONT)
    plt.xlabel('Starting Temperature', fontsize=XFONT)
    plt.ylabel('Total Abundance', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.scatter(start_temp_0, abundance_end_0, label='Dyke/Weaver Model')
    plt.scatter(start_temp_2, abundance_end_2, label="Survival Threshold 0.2")
    plt.xticks(np.arange(0,100,5))
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('abundance_alive_at_each_start_temperature_at_the_end_of_simulation.jpg')
    plt.show()

#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================


#=======================================================================================================================
data_verification()
#=======================================================================================================================
#=======================================================================================================================
#number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end()
#=======================================================================================================================
#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R()
#=======================================================================================================================
#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only()
#=======================================================================================================================
#=======================================================================================================================
#end_temperature_both_zero_zero2_with_alives_overlay()
#=======================================================================================================================
#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels()
#=======================================================================================================================
#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================