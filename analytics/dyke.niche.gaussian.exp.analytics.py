import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import statistics
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase
import pandas as pd



RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

count = 0

plt.rcParams["font.family"] = "Times New Roman"


for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.niche.gaussian.exp.data")
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

        fig, ax = plt.subplots(figsize=(20,20), dpi=200)

        plt.title('JI Model', fontsize=TFONT)
        if(each_survival_threshold == 0.2):
            plt.title('ST Model', fontsize=TFONT)
        ax.set_xlabel('Starting Temperature', fontsize=XFONT)
        ax.set_ylabel('Number of Simulations', fontsize=YFONT)
        plt.xticks(fontsize=X_TICKS)
        plt.yticks(fontsize=Y_TICKS)

        alive_below = []
        alive_above = []

        data_size = 0
        for each_start_temp in uniq_start_temps:
            index_1 = 0
            below_zero = 0
            above_zero = 0
            data_size=0
            for each_row in x:
                if(each_row == each_start_temp and z[index_1] == each_survival_threshold):
                    if(y[index_1] <= 0):
                        below_zero +=1
                    if(y[index_1] > 0):
                        above_zero +=1
                    data_size +=1
                    #print(each_survival_threshold,each_start_temp,x[index_1],y[index_1],z[index_1])
                index_1 +=1

            alive_below.append(below_zero)
            alive_above.append(above_zero)

            #main_result.append([each_start_temp,each_survival_threshold,below_zero,above_zero])

        X=[]
        for each in uniq_start_temps:
            X.append(str(each))
        plt.ylim([0, data_size + 5])

        X_axis = np.arange(len(X))

        p1 = ax.bar(X_axis, alive_above,  label='Simulations with alive species')
        p2 = ax.bar(X_axis, alive_below, bottom=alive_above , label='Simulations with no alive species')

        ax.set_xticks(X_axis)
        ax.set_xticklabels(X)
        ax.legend()

        # Label with label_type 'center' instead of the default 'edge'
        ax.bar_label(p1, label_type='center', fontsize=20)
        if(each_survival_threshold == 0.2):
            ax.bar_label(p2, label_type='center', fontsize=20)
        #ax.bar_label(p2)
        ax.legend(loc='best', fontsize=30)
        plt.tight_layout()
        plt.savefig('number_of_simulations_that_have_zero_alives_vs_more_3than_zero_alives_at_the_end_'+str(each_survival_threshold)+'.jpg')
        plt.show()



#=======================================================================================================================
#number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end()
#=======================================================================================================================
def number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end2():


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

    x = [] #start temp
    y = [] #number alive at the end
    z = [] #survival threshold

    #RESULT_DATA.append((omega,mu,niche,survival_threshold,env_start,env_end,num_alive_start,num_alive_end))

    for data_point in RESULT_DATA:
        if(data_point[2]==10):
            x.append(data_point[4])
            y.append(data_point[7])
            z.append(data_point[3])


    uniq_start_temps = np.unique(np.array(x))
    uniq_survivals = np.unique(np.array(z))
    main_result = []

    #print(uniq_start_temps)
    #print(uniq_survivals)



    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    plt.title('NW Model', fontsize=TFONT)

    ax.set_xlabel('Starting Temperature', fontsize=XFONT)
    ax.set_ylabel('Number of Simulations', fontsize=YFONT)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)

    alive_below = []
    alive_above = []

    data_size = 0
    for each_start_temp in uniq_start_temps:
        index_1 = 0
        below_zero = 0
        above_zero = 0
        data_size=0
        for each_row in x:
            if(each_row == each_start_temp):
                if(y[index_1] <= 0):
                    below_zero +=1
                if(y[index_1] > 0):
                    above_zero +=1
                data_size +=1
                #print(each_survival_threshold,each_start_temp,x[index_1],y[index_1],z[index_1])
            index_1 +=1

        alive_below.append(below_zero)
        alive_above.append(above_zero)

        #main_result.append([each_start_temp,each_survival_threshold,below_zero,above_zero])

    X=[]
    for each in uniq_start_temps:
        X.append(str(each))
    plt.ylim([0, data_size + 5])

    X_axis = np.arange(len(X))

    p1 = ax.bar(X_axis, alive_above,  label='Simulations with alive species')
    p2 = ax.bar(X_axis, alive_below, bottom=alive_above , label='Simulations with no alive species')

    ax.set_xticks(X_axis)
    ax.set_xticklabels(X)
    ax.legend()

    # Label with label_type 'center' instead of the default 'edge'
    ax.bar_label(p1, label_type='center', fontsize=20)
    ax.bar_label(p2, label_type='center', fontsize=20)
    #ax.bar_label(p2)
    ax.legend(loc='best', fontsize=30)
    plt.tight_layout()
    plt.savefig('number_of_simulations_that_have_zero_alives_vs_more_than_ze1ro_alives_at_the_end_10.jpg')
    plt.show()



#=======================================================================================================================
#number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end2()
#=======================================================================================================================


PALETTE = 'viridis'


class HandlerColormap(HandlerBase):
    def __init__(self, cmap, num_stripes=8, **kw):
        HandlerBase.__init__(self, **kw)
        self.cmap = cmap
        self.num_stripes = num_stripes
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        stripes = []
        for i in range(self.num_stripes):
            s = Rectangle([xdescent + i * width / self.num_stripes, ydescent],
                          width / self.num_stripes,
                          height,
                          fc=self.cmap((2 * i + 1) / (2 * self.num_stripes)),
                          transform=trans)
            stripes.append(s)
        return stripes


def ji_model_total_abundance_per_start_temperature():

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


    ji_start_temp = []
    ji_total_abundance_end = []

    totals = []
    UNIQ_SAMPLES = []

    for data_point in RESULT_DATA:
        if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
            UNIQ_SAMPLES.append((data_point[0],data_point[1]))

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and data_point[3]==0):
                ji_start_temp.append(data_point[4])
                ji_total_abundance_end.append(data_point[9])


    uniq_start_temps = np.unique(np.array(ji_start_temp))

##

    for uniq_sample in UNIQ_SAMPLES:
        for each_temp in uniq_start_temps:
            JI = 0
            ST = 0
            for data_point in RESULT_DATA:
                if(uniq_sample == (data_point[0],data_point[1]) and each_temp == data_point[4]):
                    if(data_point[2]==5 and data_point[3]==0):
                        if(data_point[9] < 0.1):
                            JI = 1

                    if(data_point[2]==5 and data_point[3]==0.2):
                        if(data_point[9] > 1):
                            ST = 1
            #if(JI == 1 and ST == 1):
            #    print(uniq_sample[0])
            #    print(uniq_sample[1])
            #    print(each_temp)


    zipped = list(zip(ji_start_temp, ji_total_abundance_end))
    df = pd.DataFrame(zipped, columns=['start_temp', 'abundance_end'])

    fig, ax = plt.subplots(figsize=(20,20), dpi= 200)

    #sns.boxplot(ji_start_temp, ji_total_abundance_end, palette=PALETTE)
    sns.pointplot(data = df, x = 'start_temp', y = 'abundance_end', palette=PALETTE, dodge=True, join=True, errwidth = 3, capsize = 0.5, markersize = 50)
    sns.stripplot(data = df, x = 'start_temp', y = 'abundance_end', palette=PALETTE,edgecolor='green')

    plt.title('End Abundance', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Total End Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.tight_layout()
    #cmaps = [plt.cm.spring]
    cmaps = [plt.cm.viridis]
    cmap_labels = ["JI Model"]
    # create proxy artists as handles:
    cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
    handler_map = dict(zip(cmap_handles,
                           [HandlerColormap(cm, num_stripes=8) for cm in cmaps]))
    plt.legend(handles=cmap_handles,
               labels=cmap_labels,
               handler_map=handler_map,
               fontsize=12, prop={'size': 30})

    plt.savefig('ji_model_total_abundance_per_start_temperatur9e.jpg')
    plt.show()



#=======================================================================================================================
def number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R():

    start_temp = []
    end_temp = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0)):
            start_temp.append(data_point[4])
            end_temp.append(data_point[5])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    plt.title('JI Model simulations ending with respect to the essential range', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Essential Range', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    uniq_start_temps = np.unique(np.array(start_temp))


    results_temp = []
    results_inside = []
    results_outside = []

    data_size = 0

    for each_start_temp in uniq_start_temps:
        index_1 = 0
        inside_b = 0
        outside_b = 0
        data_size = 0
        for each_row in start_temp:
            if(each_row == each_start_temp):
                if(end_temp[index_1] >= 0 and end_temp[index_1] <= 100):
                    inside_b +=1
                if(end_temp[index_1] < 0 or end_temp[index_1] > 100):
                    outside_b +=1
                if(each_row == each_start_temp):
                    data_size+=1


            index_1 +=1
        results_temp.append(each_start_temp)
        results_inside.append(inside_b)
        results_outside.append(outside_b)

    X=[]
    for each in uniq_start_temps:
        X.append(str(each))

    plt.ylim([0, data_size + 7])
    X_axis = np.arange(len(X))

    b_inside = []
    b_outside = []

    index_temp = 0

    for each_temp in X_axis:

        b_inside.append(0)
        b_outside.append(results_inside[index_temp])
        index_temp +=1


    p1 = ax.bar(X_axis, results_inside, bottom = b_inside, label = 'Simulations within the essential range')
    p2 = ax.bar(X_axis, results_outside, bottom = b_outside , label = 'Simulations outside the essential range')

    ax.set_xticks(X_axis, label = X)
    # Label with label_type 'center' instead of the default 'edge'
    ax.bar_label(p1, label_type='center', fontsize=25)
    ax.bar_label(p2, label_type='center', fontsize=25)

    #ax.bar_label(p2)
    ax.legend(loc='best', fontsize=25)
    plt.tight_layout()
    plt.savefig('number_of_simulations_that_have_end_temperature_in8side_0R_and_outside_0R_JI.jpg' )
    plt.show()

#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R()
#=======================================================================================================================

#=======================================================================================================================
def number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R_ST():

    start_temp = []
    end_temp = []
    number_alive_end = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0.2)):
            start_temp.append(data_point[4])
            end_temp.append(data_point[5])
            number_alive_end.append(data_point[7])


    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    plt.title('ST Model simulations ending with respect to the essential range', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Essential Range', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    uniq_start_temps = np.unique(np.array(start_temp))


    results_temp = []
    results_inside = []
    results_outside = []

    data_size = 0

    for each_start_temp in uniq_start_temps:
        index_1 = 0
        inside_b = 0
        outside_b = 0
        data_size = 0
        for each_row in start_temp:
            if(each_row == each_start_temp):
                if(end_temp[index_1] >= 0 and end_temp[index_1] <= 100):
                    if(number_alive_end[index_1] > 0):
                        inside_b +=1
                    else:
                        outside_b+=1
                if(end_temp[index_1] < 0 or end_temp[index_1] > 100):
                    outside_b +=1
                if(each_row == each_start_temp):
                    data_size+=1


            index_1 +=1
        results_temp.append(each_start_temp)
        results_inside.append(inside_b)
        results_outside.append(outside_b)

    X=[]
    for each in uniq_start_temps:
        X.append(str(each))

    plt.ylim([0, data_size + 7])
    X_axis = np.arange(len(X))

    b_inside = []
    b_outside = []

    index_temp = 0

    for each_temp in X_axis:

        b_inside.append(0)
        b_outside.append(results_inside[index_temp])
        index_temp +=1


    p1 = ax.bar(X_axis, results_inside, bottom = b_inside, label = 'Simulations within the essential range')
    p2 = ax.bar(X_axis, results_outside, bottom = b_outside , label = 'Simulations outside the essential range')

    ax.set_xticks(X_axis, label = X)
    # Label with label_type 'center' instead of the default 'edge'
    ax.bar_label(p1, label_type='center', fontsize=25)
    ax.bar_label(p2, label_type='center', fontsize=25)

    #ax.bar_label(p2)
    ax.legend(loc='best', fontsize=25)
    plt.tight_layout()
    plt.savefig('number_of_simulations_that_have_end_temperature_inside_0R_and_out7side_0R_ST.jpg' )
    plt.show()

#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R_ST()
#=======================================================================================================================


def function_stacked():

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
        if(data_point[2]==5 and (data_point[3]==0)):
            x.append(data_point[4])
            y.append(data_point[5])
            z.append(data_point[3])
            al.append(data_point[7])


    uniq_start_temps = np.unique(np.array(x))
    uniq_start_temps.sort()
    uniq_survivals = np.unique(np.array(z))
    uniq_survivals.sort()

    main_result = []

    #print(uniq_start_temps)
    #print(uniq_survivals)

    for each_survival_threshold in tqdm(uniq_survivals):

        fig, ax = plt.subplots(figsize=(20,20), dpi=200)

        plt.title('JI Model simulations ending with respect to the essential range', fontsize=40)
        if(each_survival_threshold == 0.2):
            plt.title('ST Model simulations ending with respect to the essential range', fontsize=40)
        ax.set_xlabel('Starting Temperature', fontsize=40)
        ax.set_ylabel('Essential Range', fontsize=40)
        plt.xticks(fontsize=X_TICKS)
        plt.yticks(fontsize=Y_TICKS)

        inside_bounds = []
        outside_bounds = []
        alives_bounds = []

        data_size = 0

        for each_start_temp in uniq_start_temps:
            index_1 = 0
            inside_b = 0
            outside_b = 0
            zero_alives = 0
            data_size=0
            for each_row in x:
                if(x[index_1] == each_start_temp and z[index_1] == each_survival_threshold and al[index_1] > 0): # al = Number of alive species greater than one
                    if(y[index_1] > 0 and y[index_1] < 100):
                        inside_b +=1
                    if(y[index_1] < 0 or y[index_1] > 100):
                        outside_b +=1
                if(x[index_1] == each_start_temp and z[index_1] == each_survival_threshold and al[index_1] == 0):
                    zero_alives +=1
                if(x[index_1] == each_start_temp and z[index_1] == each_survival_threshold):
                    data_size+=1

                index_1 +=1

            inside_bounds.append(inside_b)
            outside_bounds.append(outside_b)
            alives_bounds.append(zero_alives)

            main_result.append([each_start_temp,each_survival_threshold,inside_b,outside_b])

        X=[]
        for each in uniq_start_temps:
            X.append(str(each))

        plt.ylim([0, data_size + 7])


        X_axis = np.arange(len(X))

        #print()
        #print(inside_bounds)
        #print(outside_bounds)
        #print(alives_bounds)

        b_inside = []
        b_outside = []
        b_alive = []

        index_temp = 0

        for each_temp in X_axis:

            b_inside.append(0)
            b_outside.append(inside_bounds[index_temp])
            b_alive.append(inside_bounds[index_temp] + outside_bounds[index_temp])
            #alives_bounds[index_temp]

            index_temp +=1
        #print("====")
        #print(b_inside)
        #print(b_outside)
        #print(b_alive)

        p1 = ax.bar(X_axis, inside_bounds, bottom = b_inside, label = 'Simulations within the essential range')
        p2 = ax.bar(X_axis, outside_bounds, bottom = b_outside , label = 'Simulations outside the essential range')
        p3 = ax.bar(X_axis, alives_bounds, bottom = b_alive , label = 'Simulations with no alive species')
        ax.set_xticks(X_axis, label = X)
        # Label with label_type 'center' instead of the default 'edge'
        ax.bar_label(p1, label_type='center', fontsize=25)
        ax.bar_label(p2, label_type='center', fontsize=25)
        ax.bar_label(p3, label_type='center', fontsize=25)

        #ax.bar_label(p2)
        ax.legend(loc='best', fontsize=25)
        plt.tight_layout()
        plt.savefig('number_of_simulations_that_have_end_temperature_i5nside_0R_and_outside_0R_' + str(each_survival_threshold) + '.jpg' )
        plt.show()

def number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R2():


    start_temp = []
    end_temp = []
    number_alive_end = []

    for data_point in RESULT_DATA:
        if(data_point[2]==10):
            start_temp.append(data_point[4])
            end_temp.append(data_point[5])
            number_alive_end.append(data_point[7])


    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    plt.title('NW Model simulations ending with respect to the essential range', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Essential Range', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    uniq_start_temps = np.unique(np.array(start_temp))


    results_temp = []
    results_inside = []
    results_outside = []

    data_size = 0

    for each_start_temp in uniq_start_temps:
        index_1 = 0
        inside_b = 0
        outside_b = 0
        data_size = 0
        for each_row in start_temp:
            if(each_row == each_start_temp):
                if(end_temp[index_1] >= 0 and end_temp[index_1] <= 100):
                    if(number_alive_end[index_1] > 0):
                        inside_b +=1
                    else:
                        outside_b+=1
                if(end_temp[index_1] < 0 or end_temp[index_1] > 100):
                    outside_b +=1
                if(each_row == each_start_temp):
                    data_size+=1


            index_1 +=1
        results_temp.append(each_start_temp)
        results_inside.append(inside_b)
        results_outside.append(outside_b)

    X=[]
    for each in uniq_start_temps:
        X.append(str(each))

    plt.ylim([0, data_size + 7])
    X_axis = np.arange(len(X))

    b_inside = []
    b_outside = []

    index_temp = 0

    for each_temp in X_axis:

        b_inside.append(0)
        b_outside.append(results_inside[index_temp])
        index_temp +=1


    p1 = ax.bar(X_axis, results_inside, bottom = b_inside, label = 'Simulations within the essential range')
    p2 = ax.bar(X_axis, results_outside, bottom = b_outside , label = 'Simulations outside the essential range')

    ax.set_xticks(X_axis, label = X)
    # Label with label_type 'center' instead of the default 'edge'
    ax.bar_label(p1, label_type='center', fontsize=25)
    ax.bar_label(p2, label_type='center', fontsize=25)

    #ax.bar_label(p2)
    ax.legend(loc='best', fontsize=25)
    plt.tight_layout()
    plt.savefig('number_of_simulations_that_have_end_temperature_insid2e_0R_and_outside_0R_NW.jpg' )
    plt.show()

#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R2()
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


    fig, ax = plt.subplots(figsize=(20,20), dpi=200)
    plt.title('Essential Range', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Simulation Bounds', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)


    main_results = []

    data_size = 0

    for each_temp in tqdm(uniq_start_temps):
        both_all = []
        zero_all = []
        zero2_all = []
        data_size = 0
        zero_alives = []

        for each_sample in UNIQ_SAMPLES:
            index = 0
            both = 0
            jt_inside_0R = 0
            st_inside_0R = 0
            zeroalive = 0
            for each_data in x:
                if(x[index] == each_temp and ((omega[index], mu[index])) == each_sample and al[index] > 0):
                    if(z[index]==0 and (y[index] > 0 and y[index] < 100)):
                        jt_inside_0R +=1
                    if(z[index]==0.2 and (y[index] > 0 and y[index] < 100)):
                        st_inside_0R +=1
                if(x[index] == each_temp and ((omega[index], mu[index])) == each_sample and al[index] == 0 and z[index]==0.2):
                    zeroalive +=1

                index +=1

            if(jt_inside_0R > 0 and st_inside_0R > 0):
                both_all.append(1)
            if(jt_inside_0R > 0 and st_inside_0R == 0):
                zero_all.append(1)
            if(jt_inside_0R == 0 and st_inside_0R > 0):
                zero2_all.append(1)
            zero_alives.append(zeroalive)
            data_size+=1

        main_results.append((each_temp,sum(both_all), sum(zero_all), sum(zero2_all), sum(zero_alives)))

    bar_both = []
    bar_zero_all = []
    bar_zero2_all = []
    bar_zero_alives = []

    for each_result in main_results:
        bar_both.append(each_result[1])
        bar_zero_all.append(each_result[2])
        bar_zero2_all.append(each_result[3])
        bar_zero_alives.append(each_result[4])

    X=[]
    for each in uniq_start_temps:
        X.append(str(each))

    plt.ylim([0, data_size + 7])

    X_axis = np.arange(len(X))

    b_bar_both = []
    b_bar_zero_all = []
    b_bar_zero2_all = []
    b_bar_zero_alives = []

    index_temp = 0
    for each_temp in X_axis:
        b_bar_both.append(0)
        b_bar_zero_all.append(bar_both[index_temp])
        b_bar_zero2_all.append(bar_both[index_temp] + bar_zero_all[index_temp])
        b_bar_zero_alives.append(bar_both[index_temp] + bar_zero_all[index_temp] + bar_zero2_all[index_temp])

        index_temp +=1

    p1 = ax.bar(X_axis, bar_both, bottom = b_bar_both,  label = 'Both inside essential range')
    p2 = ax.bar(X_axis, bar_zero_all, bottom = b_bar_zero_all ,label = 'JI Model only inside range')
    p3 = ax.bar(X_axis, bar_zero2_all, bottom = b_bar_zero2_all,  label = 'ST Model only inside range')
    p4 = ax.bar(X_axis, bar_zero_alives, bottom = b_bar_zero_alives,  label = 'ST Model no alive species')

    ax.set_xticks(X_axis, label = X)
    # Label with label_type 'center' instead of the default 'edge'
    ax.bar_label(p1, label_type='center', fontsize=25)
    ax.bar_label(p2, label_type='center', fontsize=25)
    ax.bar_label(p3, label_type='center', fontsize=25)
    ax.bar_label(p4, label_type='center', fontsize=25)
    #ax.bar_label(p2)
    ax.legend(loc='best', fontsize=25)
    plt.tight_layout()
    plt.savefig('number_of_simulations_that_have_end_temperature_both_insid3e_dyke_weaver_inside_only_truncated_inside_only.jpg')
    plt.show()



#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only()
#=======================================================================================================================








###





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

    fig, ax = plt.subplots(figsize=(20,20), dpi= 200)

    #sns.boxplot(ji_start_temp, ji_total_abundance_end, palette=PALETTE)
    sns.pointplot(start_temp_2, alive_start_2, palette=PALETTE,  join=False, errwidth = 3, capsize = 0.5, markersize = 50)
    sns.stripplot(start_temp_2, alive_start_2, palette=PALETTE,edgecolor='green', dodge=True)

    plt.title('Number of alive species at the start of the simulation', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    cmaps = [plt.cm.viridis]
    cmap_labels = ["ST Model"]
    # create proxy artists as handles:
    cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
    handler_map = dict(zip(cmap_handles,
                           [HandlerColormap(cm, num_stripes=8) for cm in cmaps]))
    plt.legend(handles=cmap_handles,
               labels=cmap_labels,
               handler_map=handler_map,
               fontsize=12, prop={'size': 30})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_start_of_s1imulation.jpg')
    plt.show()


#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
def number_alive_at_each_start_temperature_at_the_start_of_simulation_mini_color_bars_working():

   #MINI COLOR BARS WORKING

    st_start_temp = []
    st_alive_start = []
    st_niche_size = []
    nw_start_temp = []
    nw_alive_start = []
    nw_niche_size = []



    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2]==5):
                st_start_temp.append(data_point[4])
                st_alive_start.append(data_point[6])
                st_niche_size.append(data_point[2])
            if(data_point[2]==10):
                nw_start_temp.append(data_point[4])
                nw_alive_start.append(data_point[6])
                nw_niche_size.append(data_point[2])

    zipped_st = list(zip(st_start_temp, st_alive_start, st_niche_size))
    df_st = pd.DataFrame(zipped_st, columns=['st_start_temp', 'st_alive_start', 'st_niche_size'])

    zipped_nw = list(zip(nw_start_temp, nw_alive_start, nw_niche_size))
    df_nw = pd.DataFrame(zipped_nw, columns=['nw_start_temp', 'nw_alive_start', 'nw_niche_size'])


    fig, ax = plt.subplots(figsize=(20,20), dpi=200)


    #colors = ['#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7', '#f7fbff']
    #customPalette = sns.set_palette(sns.color_palette(colors))
    #color1 = sns.color_palette("rocket_r", as_cmap=True)
    #color2 = sns.color_palette("viridis", as_cmap=True)
    #sns.stripplot(data = df, x = 'start_temp', y = 'alive_start', hue = 'niche_size', palette = customPalette, jitter = 0.25, dodge = True)
    #sns.pointplot(data = df, x = 'start_temp', y = 'alive_start', hue = 'niche_size', palette = customPalette, dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)

    sns.stripplot(data = df_st, x = 'st_start_temp', y = 'st_alive_start', palette = 'autumn', jitter = 0.40, dodge = True)
    sns.pointplot(data = df_st, x = 'st_start_temp', y = 'st_alive_start', palette = 'autumn', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50, join=True)
    sns.stripplot(data = df_nw, x = 'nw_start_temp', y = 'nw_alive_start', palette = 'viridis', jitter = 0.25, dodge = True)
    sns.pointplot(data = df_nw, x = 'nw_start_temp', y = 'nw_alive_start',  palette = 'viridis', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50, join=True)


    plt.legend(prop={'size': 30})
    plt.title('Number of alive species at the start of the simulation', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    cmaps = [plt.cm.autumn, plt.cm.viridis]
    cmap_labels = ["ST Model", "NW Model"]
    # create proxy artists as handles:


    cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
    handler_map = dict(zip(cmap_handles,
                           [HandlerColormap(cm, num_stripes=8) for cm in cmaps]))
    plt.legend(handles=cmap_handles,
               labels=cmap_labels,
               handler_map=handler_map,
               fontsize=12, prop={'size': 30})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_st1art_of_simulation_ST_NW.jpg')
    plt.show()




def number_alive_at_each_start_temperature_at_the_start_of_simulation2():

    start_temp = []
    alive_start = []
    niche_size = []

    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2] == 5):
                start_temp.append(data_point[4])
                alive_start.append(data_point[6])
                niche_size.append("ST Model")

            if(data_point[2] == 10):
                start_temp.append(data_point[4])
                alive_start.append(data_point[6])
                niche_size.append("NW Model")

    zipped = list(zip(start_temp, alive_start, niche_size))
    df = pd.DataFrame(zipped, columns=['start_temp', 'alive_start', 'niche_size'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    sns.stripplot(data = df, x = 'start_temp', y = 'alive_start', hue = 'niche_size', jitter = 0.25, dodge = True)
    sns.pointplot(data = df, x = 'start_temp', y = 'alive_start', hue = 'niche_size', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)
    plt.legend(prop={'size': 30})
    plt.title('Number of alive species at the start of the simulation', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)

    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], prop={'size': 25})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_start_of_s1imulation_ST_NW.jpg')
    plt.show()


#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_start_of_simulation2()
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

    start_temp_2 = []
    alive_end_2 = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                alive_end_2.append(data_point[7])


    zipped = list(zip(start_temp_2, alive_end_2))
    df = pd.DataFrame(zipped, columns=['start_temp', 'alive_end'])

    fig, ax = plt.subplots(figsize=(20,20), dpi= 200)

    #sns.boxplot(ji_start_temp, ji_total_abundance_end, palette=PALETTE)
    sns.pointplot(data = df, x = 'start_temp', y='alive_end', palette=PALETTE,  join=True, errwidth = 3, capsize = 0.5, markersize = 50)
    sns.stripplot(data = df, x = 'start_temp', y='alive_end', palette=PALETTE, edgecolor='green', dodge=True)

    plt.title('Number of alive species at the end of the simulation', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    cmaps = [plt.cm.viridis]
    cmap_labels = ["ST Model"]
    # create proxy artists as handles:
    cmap_handles = [Rectangle((0, 0), 1, 1) for _ in cmaps]
    handler_map = dict(zip(cmap_handles,
                           [HandlerColormap(cm, num_stripes=8) for cm in cmaps]))
    plt.legend(handles=cmap_handles,
               labels=cmap_labels,
               handler_map=handler_map,
               fontsize=12, prop={'size': 30})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_end_of_simul5ation.jpg')
    plt.show()



#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================
def number_alive_at_each_start_temperature_at_the_end_of_simulation2():


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



    start_temp = []
    alive_end = []
    niche_size = []

    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2] == 5):
                start_temp.append(data_point[4])
                alive_end.append(data_point[7])
                niche_size.append("ST Model")
            if(data_point[2] == 10):
                start_temp.append(data_point[4])
                alive_end.append(data_point[7])
                niche_size.append("NW Model")

    zipped = list(zip(start_temp, alive_end, niche_size))
    df = pd.DataFrame(zipped, columns=['start_temp', 'alive_end', 'niche_size'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    #color1 = 'orange'
    #color2 = 'blue'

    sns.stripplot(data = df, x = 'start_temp', y = 'alive_end', hue = 'niche_size', jitter = 0.25, dodge = True)
    sns.pointplot(data = df, x = 'start_temp', y = 'alive_end', hue = 'niche_size', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)
    plt.legend(prop={'size': 30})
    plt.title('Number of alive species at the end of the simulation', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    #red_patch = mpatches.Patch(color=color1, label='The red data')
    #red_patch = mpatches.Patch(color=color2, label='The red data')
    #plt.legend(handles=[red_patch])
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], prop={'size': 25})
    plt.tight_layout()
    plt.savefig('number_alive_at_each_start_temperature_at_the_end_of_simulatio2n_ST_NW.jpg')
    plt.show()


#=======================================================================================================================
#number_alive_at_each_start_temperature_at_the_end_of_simulation2()
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

    start_temp_2 = []
    alives_start_2 = []
    alive_end_2 = []
    difference_end_start = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0.2)):
            start_temp_2.append(data_point[4])
            alives_start_2.append(data_point[6])
            alive_end_2.append(data_point[7])
            difference_end_start.append(data_point[7] - data_point[6])


    #print(start_temp_0)
    #print(alive_end_0)

    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(start_temp_2))
    uniq_start_temps.sort()

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


    # DIFF ##########################
    temps_diff = []
    avg_diff = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += difference_end_start[index]
                count += 1

            index += 1
        temps_diff.append(each_temp)
        avg_diff.append((sum/count))



    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Average number of alive species at the start and end of the simulations', fontsize=40)
    plt.xlabel('Starting Temperature', fontsize=40)
    plt.ylabel('Average Number of Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    #plt.plot(dw_temp_s, dw_alive_s, label='Dyke/Weaver Alive Species at Start')
    #plt.plot(dw_temp_e, dw_alive_e, label='Dyke/Weaver Alive Species at End')
    plt.plot(st_temp_s, st_alive_s, label='ST Model Alive species at Start')
    plt.plot(st_temp_e, st_alive_e, label='ST Model Alive species at End')
    plt.plot(temps_diff, avg_diff, label='ST Model Average Difference')
    plt.xticks(np.arange(0,101,5))
    plt.legend(prop={'size': 25})
    plt.tight_layout()
    plt.savefig('average_number_alive_at_each_start_temperature_at_the_start_an3d_end_of_simulation.jpg')
    plt.show()


#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
def average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation2():

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

    st_start_temp = []
    st_alives_start = []
    st_alive_end = []
    st_difference_end_start = []
    nw_start_temp = []
    nw_alives_start = []
    nw_alive_end = []
    nw_difference_end_start = []

    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2]==5):
                st_start_temp.append(data_point[4])
                st_alives_start.append(data_point[6])
                st_alive_end.append(data_point[7])
                st_difference_end_start.append(data_point[7] - data_point[6])

            if(data_point[2]==10):
                nw_start_temp.append(data_point[4])
                nw_alives_start.append(data_point[6])
                nw_alive_end.append(data_point[7])
                nw_difference_end_start.append(data_point[7] - data_point[6])



    #print(start_temp_0)
    #print(alive_end_0)

    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(st_start_temp))
    uniq_start_temps.sort()

    stats_temp_2 = []
    stats_avg_alives_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in st_start_temp:
            if(each_temp == each_instance):
                sum += st_alive_end[index]
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
        for each_instance in st_start_temp:
            if(each_temp == each_instance):
                sum += st_alives_start[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_alives_2.append((sum/count))

    st_temp_s = stats_temp_2.copy()
    st_alive_s = stats_avg_alives_2.copy()

    # DIFF ##########################
    temps_diff = []
    avg_diff = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in st_start_temp:
            if(each_temp == each_instance):
                sum += st_difference_end_start[index]
                count += 1

            index += 1
        temps_diff.append(each_temp)
        avg_diff.append((sum/count))




    stats_temp_2 = []
    stats_avg_alives_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in nw_start_temp:
            if(each_temp == each_instance):
                sum += nw_alive_end[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_alives_2.append((sum/count))

    nw_temp_e = stats_temp_2.copy()
    nw_alive_e = stats_avg_alives_2.copy()

    stats_temp_2 = []
    stats_avg_alives_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in nw_start_temp:
            if(each_temp == each_instance):
                sum += nw_alives_start[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_alives_2.append((sum/count))

    nw_temp_s = stats_temp_2.copy()
    nw_alive_s = stats_avg_alives_2.copy()

    # DIFF ##########################
    nw_temps_diff = []
    nw_avg_diff = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in nw_start_temp:
            if(each_temp == each_instance):
                sum += nw_difference_end_start[index]
                count += 1

            index += 1
        nw_temps_diff.append(each_temp)
        nw_avg_diff.append((sum/count))



    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Average number of alive species at the start and end of the simulations', fontsize=40)
    plt.xlabel('Starting Temperature', fontsize=40)
    plt.ylabel('Average Number of Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.plot(st_temp_s, st_alive_s, label='ST Model Alive species at Start')
    plt.plot(st_temp_e, st_alive_e, label='ST Model Alive species at End')
    plt.plot(temps_diff, avg_diff, label='ST Model Average Difference')
    plt.plot(nw_temp_s, nw_alive_s, label='NW Model Alive species at Start')
    plt.plot(nw_temp_e, nw_alive_e, label='NW Model Alive species at End')
    plt.plot(nw_temps_diff, nw_avg_diff, label='NW Model Average Difference')
    plt.xticks(np.arange(0,101,5))
    plt.legend(prop={'size': 25})
    plt.tight_layout()
    plt.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_en2d_of_simulation2.jpg')
    plt.show()


#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation2()
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

    fig.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels_qua8d.jpg')
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
    plt.savefig('average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_level7s.jpg')
    plt.show()



#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels()
#=======================================================================================================================



def end_temperature_both_zero_zero2_with_alives_overlay():

    plt.figure(figsize=(20,20), dpi=200)
    plt.title('End temperature of JI Model vs ST Model', fontsize=40)
    plt.xlabel('ST Model End Temperature', fontsize=40)
    plt.ylabel('JI Model End Temperature', fontsize=40)
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
    cb.set_label(label="ST Model number of alive species(end)",size=YFONT)
    plt.tight_layout()
    plt.savefig('end_temperature_both_zero_zero2_with_alives_overlay6.jpg')
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

    start_temp = []
    start_abundance = []
    survival_threshold = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp.append(data_point[4])
                start_abundance.append(data_point[8])
                survival_threshold.append("JI Model")
            if(data_point[3]==0.2):
                start_temp.append(data_point[4])
                start_abundance.append(data_point[8])
                survival_threshold.append("ST Model")


    zipped = list(zip(start_temp, start_abundance, survival_threshold))
    df = pd.DataFrame(zipped, columns=['start_temp', 'start_abundance', 'survival_threshold'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    sns.stripplot(data = df, x = 'start_temp', y = 'start_abundance', hue = 'survival_threshold', jitter = 0.25, dodge = True)
    sns.pointplot(data = df, x = 'start_temp', y = 'start_abundance', hue = 'survival_threshold', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)
    plt.legend(prop={'size': 30})
    plt.title('Total Abundance of species at the start of the simulations', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], prop={'size': 25})
    plt.tight_layout()
    plt.savefig('abundance_alive_at_each_start_temperature_at_the_start_of_simulation_JI_ST5.jpg')
    plt.show()


#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
def abundance_alive_at_each_start_temperature_at_the_start_of_simulation2():


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

    start_temp = []
    start_abundance = []
    niche = []

    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2]==5):
                start_temp.append(data_point[4])
                start_abundance.append(data_point[8])
                niche.append("ST Model")
            if(data_point[2]==10):
                start_temp.append(data_point[4])
                start_abundance.append(data_point[8])
                niche.append("NW Model")


    zipped = list(zip(start_temp, start_abundance, niche))
    df = pd.DataFrame(zipped, columns=['start_temp', 'start_abundance', 'niche_size'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    sns.stripplot(data = df, x = 'start_temp', y = 'start_abundance', hue = 'niche_size', jitter = 0.25, dodge = True)
    sns.pointplot(data = df, x = 'start_temp', y = 'start_abundance', hue = 'niche_size', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)
    plt.legend(prop={'size': 30})
    plt.title('Total Abundance of species at the start of the simulations', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], prop={'size': 25})
    plt.tight_layout()
    plt.savefig('abundance_alive_at_each_start_temperature_at_the_start_of_simulation_ST_NW4.jpg')
    plt.show()


#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_start_of_simulation2()
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

    end_temp = []
    end_abundance = []
    survival_threshold = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3] == 0):
                end_temp.append(data_point[4])
                end_abundance.append(data_point[9])
                survival_threshold.append("JI Model")

            if(data_point[3] == 0.2):
                end_temp.append(data_point[4])
                end_abundance.append(data_point[9])
                survival_threshold.append("ST Model")


    zipped = list(zip(end_temp, end_abundance, survival_threshold))
    df = pd.DataFrame(zipped, columns=['end_temp', 'end_abundance', 'survival_threshold'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    sns.stripplot(data = df, x = 'end_temp', y = 'end_abundance', hue = 'survival_threshold', jitter = 0.25, dodge = True)
    sns.pointplot(data = df, x = 'end_temp', y = 'end_abundance', hue = 'survival_threshold', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)
    plt.legend(prop={'size': 30})
    plt.title('Total Abundance of species at the end of the simulations', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], prop={'size': 25})
    plt.tight_layout()
    plt.savefig('abundance_alive_at_each_start_temperature_at_the_end_of_simulation_JI_ST3.jpg')
    plt.show()

#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================
def abundance_alive_at_each_start_temperature_at_the_end_of_simulation2():


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


    end_temp = []
    end_abundance = []
    niche = []

    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2]==5):
                end_temp.append(data_point[4])
                end_abundance.append(data_point[9])
                niche.append("ST Model")

            if(data_point[2]==10):
                end_temp.append(data_point[4])
                end_abundance.append(data_point[9])
                niche.append("NW Model")


    zipped = list(zip(end_temp, end_abundance, niche))
    df = pd.DataFrame(zipped, columns=['end_temp', 'end_abundance', 'niche_size'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)

    sns.stripplot(data = df, x = 'end_temp', y = 'end_abundance', hue = 'niche_size', jitter = 0.25, dodge = True)
    sns.pointplot(data = df, x = 'end_temp', y = 'end_abundance', hue = 'niche_size', dodge = True, errwidth = 3, capsize = 0.5, markersize = 50)
    plt.legend(prop={'size': 30})
    plt.title('Total Abundance of species at the end of the simulations', fontsize=40)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[0:2], labels[0:2], prop={'size': 25})
    plt.tight_layout()
    plt.savefig('abundance_alive_at_each_start_temperature_at_the_end_of_simulation_ST_NW1.jpg')
    plt.show()

#=======================================================================================================================
#abundance_alive_at_each_start_temperature_at_the_end_of_simulation2()
#=======================================================================================================================
def average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation():

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
    abundance_start_0 = []
    abundance_end_0 = []
    start_temp_2 = []
    abundance_start_2 = []
    abundance_end_2 = []
    ji_difference_end_start = []
    st_difference_end_start = []


    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0 or data_point[3]==0.2)):
            if(data_point[3]==0):
                start_temp_0.append(data_point[4])
                abundance_start_0.append(data_point[8])
                abundance_end_0.append(data_point[9])
                ji_difference_end_start.append(data_point[9] - data_point[8])
            if(data_point[3]==0.2):
                start_temp_2.append(data_point[4])
                abundance_start_2.append(data_point[8])
                abundance_end_2.append(data_point[9])
                st_difference_end_start.append(data_point[9] - data_point[8])



    #print(start_temp_0)
    #print(alive_end_0)

    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(start_temp_0))
    uniq_start_temps.sort()

    #print(uniq_start_temps)


    # alive end stats
    stats_temp = []
    stats_avg_abundance = []

    dw_temp_s = []
    dw_abundance_s = []
    dw_temp_e = []
    dw_abundance_e = []

    st_temp_s = []
    st_abundance_s = []
    st_temp_e = []
    st_abundance_e = []


    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += abundance_end_0[index]
                count += 1

            index += 1
        stats_temp.append(each_temp)
        stats_avg_abundance.append((sum/count))

    dw_temp_e = stats_temp.copy()
    dw_abundance_e = stats_avg_abundance.copy()

    stats_temp = []
    stats_avg_abundance = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += abundance_start_0[index]
                count += 1

            index += 1
        stats_temp.append(each_temp)
        stats_avg_abundance.append((sum/count))

    dw_temp_s = stats_temp.copy()
    dw_abundance_s = stats_avg_abundance.copy()

    stats_temp_2 = []
    stats_avg_abundance_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += abundance_end_2[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_abundance_2.append((sum/count))

    st_temp_e = stats_temp_2.copy()
    st_abundance_e = stats_avg_abundance_2.copy()

    stats_temp_2 = []
    stats_avg_abundance_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += abundance_start_2[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_abundance_2.append((sum/count))

    st_temp_s = stats_temp_2.copy()
    st_abundance_s = stats_avg_abundance_2.copy()




    ########### DIFF
    ji_temps = []
    ji_diffs_avg = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += ji_difference_end_start[index]
                count += 1

            index += 1
        ji_temps.append(each_temp)
        ji_diffs_avg.append((sum/count))


    st_temps = []
    st_diffs_avg = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += st_difference_end_start[index]
                count += 1

            index += 1
        st_temps.append(each_temp)
        st_diffs_avg.append((sum/count))

    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Average Abundance at the start and end of the simulations', fontsize=40)
    plt.xlabel('Starting Temperature', fontsize=40)
    plt.ylabel('Average Abundance of Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)

    plt.plot(dw_temp_s, dw_abundance_s, label='JI Model Abundance at Start')
    plt.plot(dw_temp_e, dw_abundance_e, label='JI Model Abundance at End')
    plt.plot(st_temp_s, st_abundance_s, label='ST Model Abundance at Start')
    plt.plot(st_temp_e, st_abundance_e, label='ST Model Abundance at End')

    plt.plot(ji_temps, ji_diffs_avg, label='JI Model Average Abundance Difference')
    plt.plot(st_temps, st_diffs_avg, label='ST Model Average Abundance Difference')

    plt.xticks(np.arange(0,101,5))

    all_vals = dw_abundance_s+dw_abundance_e+st_abundance_s+st_abundance_e+ji_diffs_avg+st_diffs_avg

    plt.yticks(np.arange(math.floor(min(all_vals)), math.ceil(max(all_vals)), 1))


    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation2.jpg')
    plt.show()


#=======================================================================================================================
#average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
def average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation2():

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
    abundance_start_0 = []
    abundance_end_0 = []
    start_temp_2 = []
    abundance_start_2 = []
    abundance_end_2 = []
    ji_difference_end_start = []
    st_difference_end_start = []


    for data_point in RESULT_DATA:
        if((data_point[2]==5 and data_point[3]==0.2) or data_point[2]==10):
            if(data_point[2]==5):
                start_temp_0.append(data_point[4])
                abundance_start_0.append(data_point[8])
                abundance_end_0.append(data_point[9])
                ji_difference_end_start.append(data_point[9] - data_point[8])
            if(data_point[2]==10):
                start_temp_2.append(data_point[4])
                abundance_start_2.append(data_point[8])
                abundance_end_2.append(data_point[9])
                st_difference_end_start.append(data_point[9] - data_point[8])



    #print(start_temp_0)
    #print(alive_end_0)

    #fig.set_size_inches(3, 1.5)

    uniq_start_temps = np.unique(np.array(start_temp_0))
    uniq_start_temps.sort()

    #print(uniq_start_temps)


    # alive end stats
    stats_temp = []
    stats_avg_abundance = []

    dw_temp_s = []
    dw_abundance_s = []
    dw_temp_e = []
    dw_abundance_e = []

    st_temp_s = []
    st_abundance_s = []
    st_temp_e = []
    st_abundance_e = []


    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += abundance_end_0[index]
                count += 1

            index += 1
        stats_temp.append(each_temp)
        stats_avg_abundance.append((sum/count))

    dw_temp_e = stats_temp.copy()
    dw_abundance_e = stats_avg_abundance.copy()

    stats_temp = []
    stats_avg_abundance = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += abundance_start_0[index]
                count += 1

            index += 1
        stats_temp.append(each_temp)
        stats_avg_abundance.append((sum/count))

    dw_temp_s = stats_temp.copy()
    dw_abundance_s = stats_avg_abundance.copy()

    stats_temp_2 = []
    stats_avg_abundance_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += abundance_end_2[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_abundance_2.append((sum/count))

    st_temp_e = stats_temp_2.copy()
    st_abundance_e = stats_avg_abundance_2.copy()

    stats_temp_2 = []
    stats_avg_abundance_2 = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += abundance_start_2[index]
                count += 1

            index += 1
        stats_temp_2.append(each_temp)
        stats_avg_abundance_2.append((sum/count))

    st_temp_s = stats_temp_2.copy()
    st_abundance_s = stats_avg_abundance_2.copy()




    ########### DIFF
    ji_temps = []
    ji_diffs_avg = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_0:
            if(each_temp == each_instance):
                sum += ji_difference_end_start[index]
                count += 1

            index += 1
        ji_temps.append(each_temp)
        ji_diffs_avg.append((sum/count))


    st_temps = []
    st_diffs_avg = []

    for each_temp in uniq_start_temps:
        index = 0
        sum = 0
        count = 0
        for each_instance in start_temp_2:
            if(each_temp == each_instance):
                sum += st_difference_end_start[index]
                count += 1

            index += 1
        st_temps.append(each_temp)
        st_diffs_avg.append((sum/count))

    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Average Abundance at the start and end of the simulations', fontsize=40)
    plt.xlabel('Starting Temperature', fontsize=40)
    plt.ylabel('Average Abundance of Alive Species', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)

    plt.plot(dw_temp_s, dw_abundance_s, label='ST Model Abundance at Start')
    plt.plot(dw_temp_e, dw_abundance_e, label='ST Model Abundance at End')
    plt.plot(st_temp_s, st_abundance_s, label='NW Model Abundance at Start')
    plt.plot(st_temp_e, st_abundance_e, label='NW Model Abundance at End')

    plt.plot(ji_temps, ji_diffs_avg, label='ST Model Average Abundance Difference')
    plt.plot(st_temps, st_diffs_avg, label='NW Model Average Abundance Difference')

    plt.xticks(np.arange(0,101,5))

    all_vals = dw_abundance_s+dw_abundance_e+st_abundance_s+st_abundance_e+ji_diffs_avg+st_diffs_avg

    plt.yticks(np.arange(math.floor(min(all_vals)), math.ceil(max(all_vals)), 1))


    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.savefig('average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation21.jpg')
    plt.show()


#=======================================================================================================================
#average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation2()
#=======================================================================================================================

def end_temperature_both_zero_zero2_with_stabilization_overlay():

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

    SUB_SET = []

    for data_point in tqdm(RESULT_DATA):
        if((data_point[3] == 0 or data_point[3] == 0.2) and data_point[4] == 50 and data_point[2] == 5):
            SUB_SET.append(data_point)

    UNIQ_SAMPLES = []
    for data_point in SUB_SET:
        if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
            UNIQ_SAMPLES.append((data_point[0],data_point[1]))


    st_end_temp = []
    dw_end_temp = []

    st_in = []
    dw_in = []


    for uniq in UNIQ_SAMPLES:
        st_in_bounds = 0
        dw_in_bounds = 0
        for each in SUB_SET:
            if(((uniq[0],uniq[1])) == ((each[0],each[1]))):
                if(each[3] == 0):
                    dw_end_temp.append(each[5])
                    if(each[5] < 100 and each[5] > 0):
                        dw_in_bounds = 1
                if(each[3] == 0.2):
                    st_end_temp.append(each[5])
                    if(each[5] < 100 and each[5] > 0):
                        st_in_bounds = 1
        st_in.append(st_in_bounds)
        dw_in.append(dw_in_bounds)



    st_x_1 = []
    dw_x_1 = []

    st_x_2 = []
    dw_x_2 = []

    st_x_3 = []
    dw_x_3 = []

    index = 0
    for item in st_end_temp:
        if(st_in[index] == 1 and dw_in[index] == 1):
            #plt.scatter(st_end_temp[index], dw_end_temp[index], s=200, label="Both IN")
            st_x_1.append(st_end_temp[index])
            dw_x_1.append(dw_end_temp[index])
        index += 1

    index = 0
    for item in st_end_temp:
        if(st_in[index] == 1 and dw_in[index] == 0):
            #plt.scatter(st_end_temp[index], dw_end_temp[index], s=200, label="ST ONLY IN")
            st_x_2.append(st_end_temp[index])
            dw_x_2.append(dw_end_temp[index])

        index += 1

    index = 0
    for item in st_end_temp:
        if(st_in[index] == 0 and dw_in[index] == 1):
            #plt.scatter(st_end_temp[index], dw_end_temp[index], s=200, label="DW ONLY IN")
            st_x_3.append(st_end_temp[index])
            dw_x_3.append(dw_end_temp[index])

        index += 1
        bt = plt.scatter(st_x_1, dw_x_1, c='b', s=200, label="Both IN")
        st = plt.scatter(st_x_2, dw_x_2, c='g', s=200, label="ST ONLY IN")
        dw = plt.scatter(st_x_3, dw_x_3, c='orange', s=200, label="DW ONLY IN")

    #plt.scatter(zero_t, two_t, c=alives, cmap = 'viridis', s=200)
    #cb = plt.colorbar(orientation="vertical")
    #cb.ax.tick_params(labelsize=20)
    #cb.set_label(label="Number of alive species",size=YFONT)
    plt.legend((bt,st,dw),["Both IN" , "ST ONLY IN", "DW ONLY IN"])


    plt.tight_layout()
    plt.savefig('end_temperature_both_zero_zero2_with_stabilization_overlay2.jpg')
    plt.show()

#=======================================================================================================================
#end_temperature_both_zero_zero2_with_stabilization_overlay()
#=======================================================================================================================
def end_temperature_both_zero_zero2_with_start_temps_overlay():

    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Start and End temperature of JI Model vs ST Model', fontsize=40)
    plt.xlabel('ST Model End Temperature', fontsize=40)
    plt.ylabel('JI Model End Temperature', fontsize=40)
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
    start_temps = []

    for data_point in tqdm(RESULT_DATA):
        if((data_point[3] == 0 or data_point[3] == 0.2) and data_point[2] == 5):
            if(data_point[3] == 0):
                zero_t.append(data_point[5])
            if(data_point[3] == 0.2):
                two_t.append(data_point[5])
                start_temps.append(data_point[4])

    plt.scatter(zero_t, two_t, c=start_temps, cmap = 'plasma', s=200)
    cb = plt.colorbar(orientation="vertical")
    cb.ax.tick_params(labelsize=20)
    cb.set_label(label="Start Temperature of ST Model",size=YFONT)
    plt.tight_layout()
    plt.savefig('end_temperature_both_zero_zero2_with_start_temps_overlay1.jpg')
    plt.show()

#=======================================================================================================================
#end_temperature_both_zero_zero2_with_start_temps_overlay()
#=======================================================================================================================

def end_temperature_both_zero_zero2_with_start_temps_dw_overlay():

    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Start and End temperature of JI Model vs ST Model', fontsize=40)
    plt.xlabel('ST Model End Temperature', fontsize=40)
    plt.ylabel('JI Model End Temperature', fontsize=40)
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
    start_temps = []

    for data_point in tqdm(RESULT_DATA):
        if((data_point[3] == 0 or data_point[3] == 0.2) and data_point[2] == 5):
            if(data_point[3] == 0):
                zero_t.append(data_point[5])
                start_temps.append(data_point[4])
            if(data_point[3] == 0.2):
                two_t.append(data_point[5])


    plt.scatter(zero_t, two_t, c=start_temps, cmap = 'plasma', s=200)
    cb = plt.colorbar(orientation="vertical")
    cb.ax.tick_params(labelsize=20)
    cb.set_label(label="Start Temperature of JI Model",size=YFONT)
    plt.tight_layout()
    plt.savefig('end_temperature_both_zero_zero2_with_start_temps_dw_overlay.jpg')
    plt.show()

#=======================================================================================================================
#end_temperature_both_zero_zero2_with_start_temps_dw_overlay()
#=======================================================================================================================

def sample_space_effects_agains_optimal_growing_temperature():

    UNIQ_SAMPLES = []
    for data_point in RESULT_DATA:
        if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
            UNIQ_SAMPLES.append((data_point[0],data_point[1]))



    plt.figure(figsize=(20,20), dpi=200)
    plt.title('Sample Space effects vs optimal growing temperature', fontsize=40)
    plt.xlabel('Optimal Growing Temperature', fontsize=40)
    plt.ylabel('Effects', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)

    for sample in UNIQ_SAMPLES:
        index = 0
        for item in sample[0]:
            plt.scatter(sample[1][index],sample[0][index])
            index+=1

    plt.show()

def ji_start_end_temperature_with_abundance():


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

    start_temp = []
    end_temp = []
    end_abundance = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0)):
            start_temp.append(data_point[4])
            end_temp.append(data_point[5])
            end_abundance.append(data_point[9])

    zipped = list(zip(start_temp, end_temp, end_abundance))
    df = pd.DataFrame(zipped, columns=['start_temp', 'end_temp', 'end_abundance'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)
    ax = sns.stripplot(x="start_temp", y="end_temp", hue = 'end_abundance', data=df, palette=PALETTE, jitter=0.25)

    norm = plt.Normalize(df['end_abundance'].min(), df['end_abundance'].max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar

    plt.title('End Temperature with End Abundance', fontsize=40)
    ax.get_legend().remove()
    cbar = ax.figure.colorbar(sm)
    cbar.set_label('JI End Abundance', rotation=270, size = 40, labelpad=1)

    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('End Temperature', fontsize=40)


    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.tight_layout()
    plt.savefig('ji_start_end_temperature_with_abundance_JI.jpg')
    plt.show()


#=======================================================================================================================
#ji_start_end_temperature_with_abundance()
#=======================================================================================================================

def st_start_end_temperature_with_abundance():


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

    start_temp = []
    end_temp = []
    end_abundance = []

    for data_point in RESULT_DATA:
        if(data_point[2]==5 and (data_point[3]==0.2)):
            start_temp.append(data_point[4])
            end_temp.append(data_point[5])
            end_abundance.append(data_point[9])

    zipped = list(zip(start_temp, end_temp, end_abundance))
    df = pd.DataFrame(zipped, columns=['start_temp', 'end_temp', 'end_abundance'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)
    ax = sns.stripplot(x="start_temp", y="end_temp", hue = 'end_abundance', data=df, palette=PALETTE, jitter=0.25)

    norm = plt.Normalize(df['end_abundance'].min(), df['end_abundance'].max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar

    plt.title('End Temperature with End Abundance', fontsize=40)
    ax.get_legend().remove()
    cbar = ax.figure.colorbar(sm)
    cbar.set_label('ST End Abundance', rotation=270, size = 40, labelpad=1)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('End Temperature', fontsize=40)

    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.tight_layout()
    plt.savefig('st_start_end_temperature_with_abundance_ST.jpg')
    plt.show()


#=======================================================================================================================
#ji_start_end_temperature_with_abundance()
#=======================================================================================================================

def nw_start_end_temperature_with_abundance():


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

    start_temp = []
    end_temp = []
    end_abundance = []

    for data_point in RESULT_DATA:
        if(data_point[2]==10):
            start_temp.append(data_point[4])
            end_temp.append(data_point[5])
            end_abundance.append(data_point[9])

    zipped = list(zip(start_temp, end_temp, end_abundance))
    df = pd.DataFrame(zipped, columns=['start_temp', 'end_temp', 'end_abundance'])

    fig, ax = plt.subplots(figsize=(20,20), dpi=200)
    ax = sns.stripplot(x="start_temp", y="end_temp", hue = 'end_abundance', data=df, palette=PALETTE, jitter=0.25)

    norm = plt.Normalize(df['end_abundance'].min(), df['end_abundance'].max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])

    # Remove the legend and add a colorbar

    plt.title('End Temperature with End Abundance', fontsize=40)
    ax.get_legend().remove()
    cbar = ax.figure.colorbar(sm)
    cbar.set_label('NW End Abundance', rotation=270, size = 40, labelpad=1)
    ax.set_xlabel('Starting Temperature', fontsize=40)
    ax.set_ylabel('End Temperature', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.tight_layout()
    plt.savefig('nw_start_end_temperature_with_abundance_NW.jpg')
    plt.show()


#=======================================================================================================================
#nw_start_end_temperature_with_abundance()
#=======================================================================================================================


#=======================================================================================================================
#sample_space_effects_agains_optimal_growing_temperature()
#=======================================================================================================================
#=======================================================================================================================
#data_verification()
#=======================================================================================================================
#sample_space_effects_agains_optimal_growing_temperature()
#=======================================================================================================================
number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end()
#=======================================================================================================================
#=======================================================================================================================
number_of_simulations_that_have_zero_alives_vs_more_than_zero_alives_at_the_end2()
#=======================================================================================================================
#=======================================================================================================================
ji_model_total_abundance_per_start_temperature()
#=======================================================================================================================
#=======================================================================================================================
number_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
number_alive_at_each_start_temperature_at_the_start_of_simulation2()
#=======================================================================================================================
#=======================================================================================================================
number_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
number_alive_at_each_start_temperature_at_the_end_of_simulation2()
#=======================================================================================================================
#=======================================================================================================================
average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation2()
#=======================================================================================================================
#=======================================================================================================================
abundance_alive_at_each_start_temperature_at_the_start_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
abundance_alive_at_each_start_temperature_at_the_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
abundance_alive_at_each_start_temperature_at_the_start_of_simulation2()
#=======================================================================================================================
#=======================================================================================================================
abundance_alive_at_each_start_temperature_at_the_end_of_simulation2()
#=======================================================================================================================
#=======================================================================================================================
average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation()
#=======================================================================================================================
#=======================================================================================================================
average_number_abundance_at_each_start_temperature_at_the_start_and_end_of_simulation2()
#=======================================================================================================================
#=======================================================================================================================
number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R()
#=======================================================================================================================
number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R_ST()
#=======================================================================================================================
number_of_simulations_that_have_end_temperature_inside_0R_and_outside_0R2()
#=======================================================================================================================
#=======================================================================================================================
ji_start_end_temperature_with_abundance()
#=======================================================================================================================
#=======================================================================================================================
st_start_end_temperature_with_abundance()
#=======================================================================================================================
#=======================================================================================================================
nw_start_end_temperature_with_abundance()
#=======================================================================================================================
#=======================================================================================================================
#number_of_simulations_that_have_end_temperature_both_inside_dyke_weaver_inside_only_truncated_inside_only()
#=======================================================================================================================
#=======================================================================================================================
#average_number_alive_at_each_start_temperature_at_the_start_and_end_of_simulation_trun_levels()
#=======================================================================================================================
#=======================================================================================================================
# >>>>>end_temperature_both_zero_zero2_with_alives_overlay()
#=======================================================================================================================
#=======================================================================================================================
#end_temperature_both_zero_zero2_with_stabilization_overlay()
#=======================================================================================================================
#=======================================================================================================================
# >>>>>end_temperature_both_zero_zero2_with_start_temps_overlay()
#=======================================================================================================================
#=======================================================================================================================
# >>>>>end_temperature_both_zero_zero2_with_start_temps_dw_overlay()






