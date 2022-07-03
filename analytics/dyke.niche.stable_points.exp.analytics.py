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
from scipy import optimize


RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data_stable'
data_archives = os.listdir(data_dr)

count = 0

plt.rcParams["font.family"] = "Times New Roman"


for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.niche.stable_points.exp.data")
    #print(count)
    try:
        omega = s['omega']
        mu = s['mu']
        ji = s['JI_STABLE_POINTS']
        st = s['ST_STABLE_POINTS']
        nw = s['NW_STABLE_POINTS']

        RESULT_DATA.append((omega,mu,ji,st,nw))
        count += 1

    finally:
        s.close()

#for item in RESULT_DATA:
#    print(item[0])
#    print(item[1])
#    print(item[2])
#    print(item[3])
#    print(item[4])
#    print("------")
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

SPECIES_K   = 100                  # ----------- Number of Biotic Components
RANGE_R     = 100                  # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 1                   # ----------- Time Step3
ENV_VARS    = 1                     # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0
ENV_START=[0]

omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

########################################################################################################################
def f1(x):
    biotic_force = []
    for y in range(SPECIES_K):
        biotic_force.append(((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))) * omega[0][y])
    return(np.sum((np.array(biotic_force, dtype=float))))

def ji_stable_point_return():


    temp_range = np.arange(-25, 125, 0.01)
    biotic_values = []
    for temp in temp_range:
        biotic_values.append(f1(temp))

    plt.plot(temp_range, biotic_values)
    plt.show()

    fsolve_search_points = []
    positive_or_negative_current = 0
    # get first point ...
    first_point = 0
    index_first_point = 0
    for biotic_point in biotic_values:
        if(biotic_point>0):
            first_point = biotic_point
            break
        if(biotic_point<0):
            first_point = biotic_point
            break
        index_first_point+=1

    if(first_point > 0):
        positive_or_negative_current = 1
    if(first_point < 0):
        positive_or_negative_current = -1


    for check_index in range(index_first_point+1, len(biotic_values)):

        if(biotic_values[check_index]>0):
            if(positive_or_negative_current == -1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = +1
        if(biotic_values[check_index]<0):
            if(positive_or_negative_current == 1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = -1
    print("JI START ===================")
    print(fsolve_search_points)

    fsolve_final_roots = []
    for crossing_zero in fsolve_search_points:
        roots = optimize.fsolve(f1,crossing_zero)
        fsolve_final_roots.append(roots)

    final_stable_points = []
    for each_point in fsolve_final_roots:
        if(f1(each_point-0.01) > 0 and f1(each_point+0.01) < 0):
            final_stable_points.append(each_point)

    biotic_force = [[] for _ in range(SPECIES_K)]
    step = 0.01


    for y in range(SPECIES_K):
        for temp in temp_range:
            biotic_force[y].append((math.e) ** ((-1) * (((abs(temp-mu[0][y])) ** 2) / (2*(NICHE**2)))) * omega[0][y])

    plt.figure(figsize=(20,20), dpi=300)
    plt.title('JI', fontsize=30)
    plt.xlabel('Temperature', fontsize=30)
    plt.ylabel('Biotic Force', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)

    for _ in range(SPECIES_K):
        plt.plot(temp_range,biotic_force[_])
    for ro in final_stable_points:
        plt.axvline(ro,color='k',linewidth=5)

    plt.plot(temp_range,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4, label='Combined Biotic Force')
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.show()

    results_points = []
    for stable_p in final_stable_points:
        results_points.append(stable_p[0])

    print(set(results_points))
    print("^^^Stable Points")
    return(set(results_points))

def f1x(x):
    biotic_force = []
    truncation=0.2

    for y in range(SPECIES_K):
        aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
        if(aliveness <= truncation and aliveness >= (-1 * truncation)):
            biotic_force.append(0)
        else:
            biotic_force.append(aliveness * omega[0][y])

    return(np.sum((np.array(biotic_force, dtype=float))))

def st_stable_point_return():

    temp_range = np.arange(-25, 125, 0.01)
    biotic_values = []
    for temp in temp_range:
        biotic_values.append(f1x(temp))

    plt.plot(temp_range, biotic_values)
    plt.show()

    fsolve_search_points = []
    positive_or_negative_current = 0
    # get first point ...
    first_point = 0
    index_first_point = 0
    for biotic_point in biotic_values:
        if(biotic_point>0):
            first_point = biotic_point
            break
        if(biotic_point<0):
            first_point = biotic_point
            break
        index_first_point+=1

    if(first_point > 0):
        positive_or_negative_current = 1
    if(first_point < 0):
        positive_or_negative_current = -1


    for check_index in range(index_first_point+1, len(biotic_values)):

        if(biotic_values[check_index]>0):
            if(positive_or_negative_current == -1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = +1
        if(biotic_values[check_index]<0):
            if(positive_or_negative_current == 1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = -1


    print("ST START ===================")
    print(fsolve_search_points)

    fsolve_final_roots = []
    for crossing_zero in fsolve_search_points:
        roots = optimize.fsolve(f1x,crossing_zero)
        fsolve_final_roots.append(roots)

    final_stable_points = []
    for each_point in fsolve_final_roots:
        if(f1x(each_point-0.01) > 0 and f1x(each_point+0.01) < 0):
            final_stable_points.append(each_point)

    temperatures = []
    alive_value = [[] for _ in range(SPECIES_K)]
    step = 0.01

    truncation = 0.2
    for y in range(SPECIES_K):
        for x in temp_range:
            aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
            if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                alive_value[y].append(0)
            else:
                alive_value[y].append(aliveness * omega[0][y])

    plt.figure(figsize=(20,20), dpi=300)
    plt.title('ST', fontsize=30)
    plt.xlabel('Temperature', fontsize=30)
    plt.ylabel('Biotic Force', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    for _ in range(SPECIES_K):
        plt.plot(temp_range,alive_value[_])

    for ro in final_stable_points:
        plt.axvline(ro,color='k',linewidth=5)

    plt.plot(temp_range,np.sum((np.array(alive_value, dtype=float)), axis=0), lw=4, label='Combined Biotic Force')
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.show()

    results_points = []
    for stable_p in final_stable_points:
        results_points.append(stable_p[0])
    print(set(results_points))
    print("^^^Stable Points")
    return(set(results_points))

########################################################################################################################

def fYa(Xe, Ni, u):
    return (((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2))))))

def fXe(Ya, Ni, u):
    return (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))) + u

def fYaI(Xe, Ni, u, T):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if(abundance <= T):
        abundance = 0

    return(abundance)

def fXe(Ya, Ni, u):
    return (u + (math.sqrt(((-1 * math.log(Ya,math.e)) * (2*(Ni**2))))))

def fXe_negative(Ya, Ni, u):
    return (u - (math.sqrt(((-1 * math.log(Ya,math.e)) * (2*(Ni**2))))))

def fYaIx(Xe, Ni, u, NRange):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if((Xe >= u + NRange) or (Xe <= u - NRange)):
        abundance = 0

    return(abundance)
########################################################################################################################

def f1x2(x):

    biotic_force = []
    truncation=0.2

    for y in range(SPECIES_K):
        NRange = (fXe(0.2, 5, mu[0][y]) - mu[0][y])
        aliveness = fYaIx(x, 10, mu[0][y], NRange)
        biotic_force.append(aliveness * omega[0][y])

    return(np.sum((np.array(biotic_force, dtype=float))))

def nw_stable_point_return():


    temp_range = np.arange(-25, 125, 0.01)
    biotic_values = []
    for temp in temp_range:
        biotic_values.append(f1x2(temp))

    plt.plot(temp_range, biotic_values)
    plt.show()

    fsolve_search_points = []
    positive_or_negative_current = 0
    # get first point ...
    first_point = 0
    index_first_point = 0
    for biotic_point in biotic_values:
        if(biotic_point>0):
            first_point = biotic_point
            break
        if(biotic_point<0):
            first_point = biotic_point
            break
        index_first_point+=1

    if(first_point > 0):
        positive_or_negative_current = 1
    if(first_point < 0):
        positive_or_negative_current = -1


    for check_index in range(index_first_point+1, len(biotic_values)):

        if(biotic_values[check_index]>0):
            if(positive_or_negative_current == -1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = +1
        if(biotic_values[check_index]<0):
            if(positive_or_negative_current == 1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = -1


    print("NW START ===================")
    print(fsolve_search_points)

    print("Reducing First ...")

    reduced_points = []
    for each_point in fsolve_search_points:
        if(f1x2(each_point-0.01) > 0 and f1x2(each_point+0.01) < 0):
            reduced_points.append(each_point)
    print("Reduced Points : " + str(reduced_points))

    fsolve_final_roots = []
    for crossing_zero in reduced_points:
        roots = optimize.fsolve(f1x2,crossing_zero)
        fsolve_final_roots.append(roots)
    print("FSolved Roots : " + str(fsolve_final_roots))

    final_stable_points = reduced_points

    temperatures = []
    alive_value = [[] for _ in range(SPECIES_K)]
    step = 0.01

    truncation = 0.2
    for x in np.arange (-25, RANGE_R+25, step):
        temperatures.append(x)

    for y in range(SPECIES_K):
        for x in np.arange (-25, RANGE_R+25, step):
            NRange = (fXe(0.2, 5, mu[0][y]) - mu[0][y])
            aliveness = fYaIx(x, 10, mu[0][y], NRange)
            alive_value[y].append(aliveness * omega[0][y])

    plt.figure(figsize=(20,20), dpi=300)
    plt.title('NW Final', fontsize=30)
    plt.xlabel('Temperature', fontsize=30)
    plt.ylabel('Biotic Force', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    #for _ in range(SPECIES_K):
    #    plt.plot(temperatures,alive_value[_])

    for ro in final_stable_points:
        plt.axvline(ro,color='k',linewidth=3)

    plt.plot(temperatures,np.sum((np.array(alive_value, dtype=float)), axis=0), lw=1, label='Combined Biotic Force')
    plt.legend(prop={'size': 30})

    x = np.arange (-25, RANGE_R+25, 5)
    plt.xticks(x)
    plt.tight_layout()
    plt.show()

    print(set(final_stable_points))
    print("^^^Stable Points")
    return(set(final_stable_points))


def all_stable_points():

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


    for uniq_s in RESULT_DATA:

        global omega
        omega = uniq_s[0]
        global mu
        mu = uniq_s[1]

        ji = ji_stable_point_return()
        st = st_stable_point_return()
        nw = nw_stable_point_return()

        print("JI: "+str(uniq_s[2]))
        print("ST: "+str(uniq_s[3]))
        print("NW: "+str(uniq_s[4]))

#all_stable_points()

def stable_points_for_ji_st_nw():

    stable_point_temp = []
    model_name = []

    for entry in RESULT_DATA:
        for sp in entry[2]:
            stable_point_temp.append(sp)
            model_name.append("JI Stable Points")
        for sp in entry[3]:
            stable_point_temp.append(sp)
            model_name.append("ST Stable Points")
        for sp in entry[4]:
            stable_point_temp.append(sp)
            model_name.append("NW Stable Points")


    zipped = list(zip(model_name, stable_point_temp))
    df = pd.DataFrame(zipped, columns=['model_name', 'stable_point_temp'])

    fig, ax = plt.subplots(figsize=(20,20), dpi= 200)

    #sns.pointplot(data = df, x = 'start_temp', y = 'abundance_end', dodge=True, join=True, errwidth = 3, capsize = 0.5, markersize = 50)
    sns.stripplot(data = df, x = 'model_name', y = 'stable_point_temp', palette = 'viridis',  edgecolor='green', dodge=True)

    plt.title('Stable Point', fontsize=40)
    ax.set_xlabel('JI, ST, NW', fontsize=40)
    ax.set_ylabel('Stable Points', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.tight_layout()
    plt.savefig('stable_points_for_ji_st_nw.jpg')
    plt.show()

stable_points_for_ji_st_nw()

def stable_point_counts_for_ji_st_nw():

    average_stable_points = []
    model_name = []

    #checks
    ji_sp = []
    st_sp = []
    nw_sp = []


    i = 0
    index = []



    for entry in RESULT_DATA:

        average_stable_points.append(len(entry[2]))
        ji_sp.append(len(entry[2]))
        model_name.append("JI Stable Points")

        average_stable_points.append(len(entry[3]))
        st_sp.append(len(entry[3]))
        model_name.append("ST Stable Points")

        average_stable_points.append(len(entry[4]))
        nw_sp.append(len(entry[4]))
        model_name.append("NW Stable Points")



    zipped = list(zip(model_name, average_stable_points))
    df = pd.DataFrame(zipped, columns=['model_name', 'average_stable_points'])

    fig, ax = plt.subplots(figsize=(20,20), dpi= 200)


    sns.stripplot(data = df, x = 'model_name', y = 'average_stable_points', palette = 'viridis',  edgecolor='green')
    sns.pointplot(data = df, x = 'model_name', y = 'average_stable_points', palette = 'autumn', errwidth = 3, capsize = 0.5, markersize = 50)

    print("JI Average Stable Points : " + str(statistics.mean(ji_sp)))
    print("ST Average Stable Points : " + str(statistics.mean(st_sp)))
    print("NW Average Stable Points : " + str(statistics.mean(nw_sp)))



    plt.title('Stable Point', fontsize=40)
    ax.set_xlabel('JI, ST, NW', fontsize=40)
    ax.set_ylabel('Stable Points', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    plt.tight_layout()
    plt.savefig('stable_point_counts_for_ji_st_nw.jpg')
    plt.show()

stable_point_counts_for_ji_st_nw()

