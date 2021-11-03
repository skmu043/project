import sys
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D


if int(len(sys.argv)) != int(2):
    print("Args: shelve file name which contains all of > (K, R, P, E, start, end, step, EN, OE, LP_Z, RUN_ID)")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN=2, OE=5, LP_Z = (10 - 100), RUN_ID : epoch")
    print("exit")
    sys.exit()

s = shelve.open(str(sys.argv[1]))

try:
    SAMPLE_SIZE = s['SAMPLE_SIZE']
    SAMPLE_STEP = s['SAMPLE_STEP']
    RUN_ID = s['RUN_ID']

    biotic_components_K = s['biotic_components_K']
    essential_range_R = s['essential_range_R']
    external_perturbation_rate_P = s['external_perturbation_rate_P']
    time_start = s['time_start']
    time_end = s['time_end']
    time_step = s['time_step']
    environment_components_N = s['environment_components_N']
    truncated_gaussian_ROUND = s['truncated_gaussian_ROUND']
    niche_width = s['niche_width']
    local_population_size = s['local_population_size']
    affects_w = s['affects_w']
    optimum_condition_u = s['optimum_condition_u']
    biotic_force_F = s['biotic_force_F']
    local_population_index = s['local_population_index']

    #global_start_temp = s['global_start_temp']
    #local_start_temp = s['local_start_temp']

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']

finally:
    s.close()

S_STEP = SAMPLE_STEP
K = biotic_components_K
R = essential_range_R
P = external_perturbation_rate_P
start = time_start
end = time_end
step = time_step
N = environment_components_N
#E = [global_start_temp, local_start_temp]
F = biotic_force_F


ROUND = truncated_gaussian_ROUND

OEn = niche_width
OE = [OEn for _ in range(K)]

w = affects_w
u = optimum_condition_u

rAx = [[] for x in range(K)]
rAxR = [[] for x in range(K)]
rNumberAlive = [[] for x in range(K)]
alpha = [[] for _ in range(K)]

def soon():
    for _ in range(K):
        al = []
        for ai in range(N):
            al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

            new_alpha = 0
            if _ in local_population_index:
                new_alpha = np.prod(al) # LOCAL
            else:
                new_alpha = al[0]       # GLOBAL

            alpha[_].append(new_alpha)

            rAx[_].append(new_alpha)
            rAxR[_].append(new_alpha * R)

            if(new_alpha > 0):
                rNumberAlive[_].append(1)
            else:
                rNumberAlive[_].append(0)

def biotic_effect_global_population():

    temperatures=[]
    for temp_value in np.arange (-50, R+50, step):
        temperatures.append(temp_value)

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force at Temperature : Global Species', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    super_biotic_force = []
    biotic_force = []

    for each_species in range(K):
        if each_species not in local_population_index:
            biotic_force_run = []
            for temp_value in np.arange(-50,R+50,step):
                biotic_force_run.append( round((math.e) ** ((-1) * (((abs((temp_value)-u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND) * w[0][each_species])
            biotic_force.append(biotic_force_run)

    for global_species in biotic_force:
        plt.plot(temperatures,global_species)

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
    super_biotic_force.append(np.sum((np.array(biotic_force, dtype=float)), axis=0))

    plt.savefig("biotic_effect_global_population.png", format = "png", dpi = 100)
    #plt.show()

mpl.use('agg') #for the 3D magic
##############################################################################################################################

def biotic_effect_global_population_extended_local():

    temperatures=[]
    for temp_value in np.arange (-50, R+50, step):
        temperatures.append(temp_value)

    temperatures_local_e=[]
    for temp_value in np.arange (-50, R+50, step):
        temperatures_local_e.append(temp_value)


    fig = plt.figure(figsize=(20,20), dpi=500)
    ax = fig.gca(projection='3d')
    ax.set_title(label = "Biotic Force at Temperature : Global Species with El")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Species Biotic Force')
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])

    super_biotic_force = []
    biotic_force = []

    for each_species in range(K):
        if each_species not in local_population_index:
            biotic_force_run = []
            for temp_value in np.arange(-50,R+50,step):
                biotic_force_run.append( round((math.e) ** ((-1) * (((abs((temp_value)-u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND) * w[0][each_species])
            biotic_force.append(biotic_force_run)



    global_sum_temp_biotic = [[],[],[]] # local temp + global temp + biotic force
    biotic_index = 0
    for each_global_temp_unit in temperatures:

        for each_local_temp_unit in temperatures_local_e:
            sum_biotic_force = 0
            for curve in biotic_force:
                sum_biotic_force += curve[biotic_index]

            global_sum_temp_biotic[0].append(each_local_temp_unit)
            global_sum_temp_biotic[1].append(each_global_temp_unit)
            global_sum_temp_biotic[2].append(sum_biotic_force)

            ax.scatter(each_local_temp_unit,each_global_temp_unit,sum_biotic_force, s=1, color="aquamarine")
        biotic_index+=1


    super_biotic_force.append(np.sum((np.array(biotic_force, dtype=float)), axis=0))

    plt.savefig("biotic_effect_global_population_local_e.png", format = "png", dpi = 100)
    #plt.show()

##############################################################################################################################



def biotic_effect_local_and_globalsum():

    # Heat Map

    fig = plt.figure(figsize=(10,10), dpi=500)
    ax = fig.gca(projection='3d')
    ax.set_title(label = "Total Biotic Force at EL + EG (Global Local Combined)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Total Biotic Force')
    ax.set_xlim([-50,R+50])
    ax.set_ylim([-50,R+50])


    heatmap = [[0 for _ in np.arange(-50,R+50,step)] for _ in np.arange(-50,R+50,step)]

    ########### LOCAL

    for each_species in range(K):
        if each_species in local_population_index:
            for global_var_temp in np.arange(-50,R+50,step):      # Env Var that affects Global and Local
                for local_var_temp in np.arange(-50,R+50,step):  # Env Var that affects Local Only
                    #print("BioticRun : ", global_var_temp,":",local_var_temp)
                    biotic_effect_local = (
                            (round((math.e) ** ((-1) * (((abs((global_var_temp)-u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND))
                            *
                            (round((math.e) ** ((-1) * (((abs((local_var_temp)-u[1][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND))
                            *
                            w[1][each_species]
                    )

                    heatmap[int(local_var_temp)][int(global_var_temp)] += biotic_effect_local


    ##### GLOBAL
    temperatures=[]
    for temp_value in np.arange (-50, R+50, step):
        temperatures.append(temp_value)

    temperatures_local_e=[]
    for temp_value in np.arange (-50, R+50, step):
        temperatures_local_e.append(temp_value)

    biotic_force = []
    for each_species in range(K):
        if each_species not in local_population_index:
            biotic_force_run = []
            for temp_value in np.arange(-50,R+50,step):
                biotic_force_run.append( round((math.e) ** ((-1) * (((abs((temp_value)-u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND) * w[0][each_species])
            biotic_force.append(biotic_force_run)

    biotic_index = 0
    for each_global_temp_unit in temperatures:
        print(biotic_index)
        for each_local_temp_unit in temperatures_local_e:
            sum_biotic_force = 0
            for curve in biotic_force:
                sum_biotic_force += curve[biotic_index]
            heatmap[int(each_local_temp_unit)][int(each_global_temp_unit)] += sum_biotic_force

        biotic_index+=1

    number_of_points = 0

    biotic_force_run = []
    global_temperatures = []
    local_temperatures = []

    for global_var_temp in np.arange(-50,R+50,step):      # Env Var that affects Global and Local
        for local_var_temp in np.arange(-50,R+50,step):  # Env Var that affects Local Only
            #print("PlotRun : ", global_var_temp,":",local_var_temp)
            if(heatmap[int(local_var_temp)][int(global_var_temp)] > 0.01 or heatmap[int(local_var_temp)][int(global_var_temp)] < -0.01):
                #ax.scatter(local_var_temp, global_var_temp, heatmap[int(local_var_temp)][int(global_var_temp)], s=1, color="yellowgreen")
                #print("Points : ", number_of_points)
                number_of_points += 1

                global_temperatures.append(global_var_temp)
                local_temperatures.append(local_var_temp)
                biotic_force_run.append(heatmap[int(local_var_temp)][int(global_var_temp)])

    ax.plot_trisurf(local_temperatures, global_temperatures, biotic_force_run,linewidth=0.2, antialiased=True,cmap=mpl.cm.jet)

    plt.savefig("biotic_effect_local_and_globalsum.png", format = "png", bbox_inches="tight")
    ax.view_init(elev=10., azim=100)
    plt.savefig("biotic_effect_local_and_globalsum_r.png", format = "png", bbox_inches="tight")


    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Global + Local Species Biotic Force', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-50, R+50)
    plt.xlim(-50, R+50)


    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('biotic_effect_local_and_globalsum_heatmap.png')

#plt.show()



def biotic_effect_local_population():

    fig = plt.figure(figsize=(10,10), dpi=500)
    ax = fig.gca(projection='3d')
    ax.set_title(label = "Biotic Force at Temperature : Local Species (affected by Eg and El)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Species Biotic Force')
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])


    for each_species in range(K):
        if each_species in local_population_index:
            biotic_force_run = []
            global_temperatures = []
            local_temperatures = []
            for global_var_temp in np.arange(-50,R+50,step):      # Env Var that affects Global and Local
                for local_var_temp in np.arange(-50,R+50,step):  # Env Var that affects Local Only

                    biotic_effect_local = (
                            (round((math.e) ** ((-1) * (((abs((global_var_temp)-u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND))
                                *
                            (round((math.e) ** ((-1) * (((abs((local_var_temp)-u[1][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND))
                                *
                            w[1][each_species]
                    )

                    if(biotic_effect_local > 0.01 or biotic_effect_local < -0.01):
                        global_temperatures.append(global_var_temp)
                        local_temperatures.append(local_var_temp)
                        biotic_force_run.append (biotic_effect_local)

            ax.scatter(local_temperatures, global_temperatures, biotic_force_run, s=1)

    plt.savefig("biotic_effect_local_population.png", format = "png" , bbox_inches="tight")
    ax.view_init(elev=10., azim=100)
    plt.savefig("biotic_effect_local_population_r.png", format = "png" , bbox_inches="tight")

    #plt.show()


    #plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
    #super_biotic_force.append(np.sum((np.array(biotic_force, dtype=float)), axis=0))



def biotic_effect_local_sum():

    # Heat Map

    fig = plt.figure(figsize=(10,10), dpi=500)
    ax = fig.gca(projection='3d')
    ax.set_title(label = "Total Biotic Force at Temperature : Local Species (affected by Eg and El)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Total Biotic Force')
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])


    heatmap = [[0 for _ in np.arange(-50,R+50,step)] for _ in np.arange(-50,R+50,step)]

    for each_species in range(K):
        if each_species in local_population_index:
            for global_var_temp in np.arange(-50,R+50,step):      # Env Var that affects Global and Local
                for local_var_temp in np.arange(-50,R+50,step):  # Env Var that affects Local Only
                    #print("BioticRun : ", global_var_temp,":",local_var_temp)
                    biotic_effect_local = (
                            (round((math.e) ** ((-1) * (((abs((global_var_temp)-u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND))
                            *
                            (round((math.e) ** ((-1) * (((abs((local_var_temp)-u[1][each_species])) ** 2) / (2*(OE[each_species]**2)))),ROUND))
                            *
                            w[1][each_species]
                    )

                    heatmap[int(local_var_temp)][int(global_var_temp)] += biotic_effect_local


    number_of_points = 0

    biotic_force_run = []
    global_temperatures = []
    local_temperatures = []

    for global_var_temp in np.arange(-50,R+50,step):      # Env Var that affects Global and Local
        for local_var_temp in np.arange(-50,R+50,step):  # Env Var that affects Local Only
            #print("PlotRun : ", global_var_temp,":",local_var_temp)
            if(heatmap[int(local_var_temp)][int(global_var_temp)] > 0.01 or heatmap[int(local_var_temp)][int(global_var_temp)] < -0.01):
                #ax.scatter(local_var_temp, global_var_temp, heatmap[int(local_var_temp)][int(global_var_temp)], s=1, color="yellowgreen")
                #print("Points : ", number_of_points)
                number_of_points += 1

                global_temperatures.append(global_var_temp)
                local_temperatures.append(local_var_temp)
                biotic_force_run.append(heatmap[int(local_var_temp)][int(global_var_temp)])

    ax.plot_trisurf(local_temperatures, global_temperatures, biotic_force_run,linewidth=0.2, antialiased=True,cmap=mpl.cm.jet)

    plt.savefig("biotic_effect_local_sum.png", format = "png", bbox_inches="tight")
    ax.view_init(elev=10., azim=100)
    plt.savefig("biotic_effect_local_sum_r.png", format = "png", bbox_inches="tight")


    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Local Species Biotic Force', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)


    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('local_species_biotic_force_heatmap.png')

#plt.show()




biotic_effect_global_population()

biotic_effect_local_population()

biotic_effect_local_sum()

biotic_effect_global_population_extended_local()

biotic_effect_local_and_globalsum()