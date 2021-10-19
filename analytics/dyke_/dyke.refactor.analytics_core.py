import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from multiprocessing import Process, Pool
import gc
from itertools import product

E_prime             =[]
F_prime             =[]
alpha_prime         =[]
rF_prime            =[]
rE_prime            =[]
rAx_prime           =[]
rAxR_prime          =[]
rAxS_prime          =[]
time_prime          =[]
biotic_force_prime  =[]
temperatures_prime  =[]
#simulation_run      =[]
rNumberAlive_prime  =[]

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)
total_to_be_processed = len(data_archives)

biotic_components_K = 0
essential_range_R = 0
time_end = 0
time_step = 0
environment_components_N = 0
truncated_gaussian_ROUND = 0
niche_width = 0
local_population_size = 0
affects_w = []
optimum_condition_u = []
biotic_force_F  = 0
exp_name = ""
data_directory = ""
shelve_file = ""
time_step  = 0
local_population_index = 0
global_start_temp = 0
local_start_temp = 0
print(gc.isenabled())
RUN_ID=0
rF = []
rE = []
time = []
OE = []
rAx = []
rAxR = []
rNumberAlive = []

def get_shelve(file):
    print("Processing  : ",print(len(rAxR_prime)), " >>> ", file)
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor_core.data")
    return(s)

def read_files_parallel(file):
    print("Processing  : ",print(len(rAxR_prime)), " >>> ", file)
    gc.collect()
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor_core.data")
    global biotic_components_K , essential_range_R, time_end, time_step, environment_components_N, \
        truncated_gaussian_ROUND, niche_width, local_population_size, affects_w, optimum_condition_u, \
        biotic_force_F, exp_name, data_directory, shelve_file, time_step, local_population_index, \
        global_start_temp, local_start_temp, RUN_ID, rF, rE, time, OE, rAx, rAxR, rNumberAlive

    try:
        #SAMPLE_SIZE = s['SAMPLE_SIZE']
        #SAMPLE_STEP = s['SAMPLE_STEP']
        RUN_ID = s['RUN_ID']

        biotic_components_K = s['biotic_components_K']
        essential_range_R = s['essential_range_R']
        #external_perturbation_rate_P = s['external_perturbation_rate_P']
        #time_start = s['time_start']
        time_end = s['time_end']
        time_step = s['time_step']
        environment_components_N = s['environment_components_N']
        truncated_gaussian_ROUND = s['truncated_gaussian_ROUND']
        niche_width = s['niche_width']
        local_population_size = s['local_population_size']
        affects_w = s['affects_w']
        optimum_condition_u = s['optimum_condition_u']
        biotic_force_F = s['biotic_force_F']

        exp_name = s['exp_name']
        data_directory = s['data_directory']
        shelve_file = s['shelve_file']

        time_step = s['time_step']
        local_population_index = s['local_population_index']

        global_start_temp = s['global_start_temp']
        local_start_temp = s['local_start_temp']


        rAx = s['rAx']
        rAxR = s['rAxR']
        rNumberAlive = s['rNumberAlive']
        alpha = s['alpha']
        rF = s['rF']
        rE = s['rE']
        time = s['time']

        OE = s['OE']

        #E_prime.append([global_start_temp, local_start_temp])
        #F_prime.append(F)
        for _ in range(biotic_components_K):
            del alpha[_][-1]
        for _ in range(biotic_components_K):
            del rAx[_][-1]
        for _ in range(biotic_components_K):
            del rAxR[_][-1]
        for _ in range(biotic_components_K):
            del rNumberAlive[_][-1]


        #alpha_prime.append(alpha)
        #rF_prime.append(rF)
        #rE_prime.append(rE)
        #rNumberAlive_prime.append(rNumberAlive)
        #time_prime.append(time)
        #rAx_prime.append(rAx)
        #rAxR_prime.append(rAxR)

    finally:
        s.close()

w = affects_w
u = optimum_condition_u

#S_STEP = SAMPLE_STEP
K = biotic_components_K

R = essential_range_R
print("R", R)
#P = external_perturbation_rate_P
#start = time_start
#end = time_end
step = time_step
N = environment_components_N
E = [0, 0]
F = biotic_force_F

ROUND = truncated_gaussian_ROUND

OEn = niche_width


def stable_points_space_heat():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Heat Map of Global and Local Species Distribution', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    # Species surviving
    index_global=0

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for optimum_condition in u[0]: # Global
        plot_x=[]
        plot_y=[]
        if(index_global not in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1

    index_global=0
    for optimum_condition in u[1]: # Local
        plot_x=[]
        plot_y=[]
        if(index_global in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')

    plt.show()

#stable_points_space_heat()

temperatures = []
def plot_alphas():

    for x in np.arange (-50, R+50, step):
        temperatures.append(x)

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force over Temperature', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    super_biotic_force = []

    for each_env_var in range(N):
        biotic_force = [[] for _ in range(K)]
        for y in range(K):
            for x in np.arange (-50, R+50, step):
                biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[each_env_var][y])) ** 2) / (2*(OE[y]**2)))) * w[each_env_var][y])


        for _ in range(K):
            plt.plot(temperatures,biotic_force[_])


        plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
        super_biotic_force.append(np.sum((np.array(biotic_force, dtype=float)), axis=0))

    sum = []
    #for _ in range(time_end/time_step):
    #    sum.append(super_biotic_force[0][_] + super_biotic_force[1][_])
    #plt.plot(temperatures, sum, lw=10)
    #plt.show()
#plot_alphas()


def stable_points_space_heat_global():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Heat Map of Global Species Distribution', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    # Species surviving
    index_global=0

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for optimum_condition in u[0]: # Global
        plot_x=[]
        plot_y=[]
        if(index_global not in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')

    plt.show()

#stable_points_space_heat_global()

def stable_points_space_heat_local():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Heat Map of Local Species Distribution', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    # Species surviving
    index_global=0

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for optimum_condition in u[1]: # Local
        plot_x=[]
        plot_y=[]
        if(index_global in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')

    plt.show()

#stable_points_space_heat_local()

def stable_points_space_final_abundance_rgb_heat():

    R = 100
    biotic_components_K = 100

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Final Abundance HeatMap', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    heatmap = [[0 for _ in range(R+10)] for _ in range(R+10)]

    index = 1
    for file in data_archives:

        s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor_core.data")

        print("Processing: ", index, "/", len(data_archives)," >>> " , file)
        index +=1
        try :
            row = s['rE']
            decompose = s['rAx']

            for _ in range(biotic_components_K):
                del decompose[_][-1]

            c_r = int(row[1][-1])
            c_g = int(row[0][-1])

            row_abundance = []

            run_length = len(decompose[0])

            for x in range(run_length):
                current_sum = 0
                for item in decompose:
                    current_sum += item[x]
                row_abundance.append(current_sum)

            if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
                #plt.scatter(row[1][0], row[0][0],s=50, marker='.', color=(float(row_abundance[0]/100), float(1), float(1)))

                heatmap[int(row[1][0])][int(row[0][0])] = row_abundance[-1]
                #plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
            else:
                #plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
                #plt.scatter(row[1][0], row[0][0], row_abundance[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
                #print(row_abundance[-1]/100)
                #plt.scatter(row[1][0], row[0][0], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(row_abundance[-1]/100)))
                #plt.scatter(row[1][0], row[0][0], marker='*', s=50, color=(float(0), float(row_abundance[-1]/100), float(0)))

                heatmap[int(row[1][0])][int(row[0][0])] = row_abundance[-1]
        finally:
            s.close()

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')
    #plt.imsave("abundanceheat.png",heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('abundanceheat.png')
    plt.show()

#stable_points_space_final_abundance_rgb_heat()

def stable_points_space_final_alive_rgb_heat():

    stable_locations = []
    R = 100
    biotic_components_K = 100
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Final Alive HeatMap', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    index = 1
    for file in data_archives:

        s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor_core.data")
        print("Processing: ", index, "/", len(data_archives)," >>> " , file)
        index +=1
        try:
            row = s['rE']
            decompose = s['rNumberAlive']

            for _ in range(biotic_components_K):
                del decompose[_][-1]

            c_r = int(row[1][-1])
            c_g = int(row[0][-1])

            row_alive = []

            run_length = len(decompose[0])

            for x in range(run_length):
                current_sum = 0
                for item in decompose:
                    current_sum += item[x]
                row_alive.append(current_sum)


            if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
                #plt.scatter(row[1][0], row[0][0],s=50, marker='.', color=(float(row_abundance[0]/100), float(1), float(1)))
                heatmap[int(row[1][0])][int(row[0][0])] = row_alive[-1]
                #plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
            else:
                #plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
                #plt.scatter(row[1][0], row[0][0], row_abundance[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
                #print(row_abundance[-1]/100)
                #plt.scatter(row[1][0], row[0][0], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(row_abundance[-1]/100)))
                #plt.scatter(row[1][0], row[0][0], marker='*', s=50, color=(float(0), float(row_abundance[-1]/100), float(0)))
                heatmap[int(row[1][0])][int(row[0][0])] = row_alive[-1]

        finally:
            s.close()

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')
    #plt.imsave("aliveheat.png",heatmap, cmap='hot', interpolation='nearest', origin='lower')

    plt.savefig('aliveheat.png')
    plt.show()

#stable_points_space_final_alive_rgb_heat()

def stable_points_global_local():

    R = 100

    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Regions', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    index = 1

    for file in data_archives:

        print("Processing: ", index, "/", len(data_archives)," >>> " , file)
        index+=1

        s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor_core.data")

        try:
            row = s['rE']

            stable_locations = []

            if((int(row[1][-1]),int(row[0][-1])) not in stable_locations):
                stable_locations.append((int(row[1][-1]),int(row[0][-1])))

            c_r = int(row[1][-1])
            c_g = int(row[0][-1])

            if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
                # float color must be between 0 and 1
                # trajectories outside 0 R
                plt.plot(row[1][0], row[0][0], marker='.', markersize = "10", color=(float(0), float(0), float(1))) # Plots Start but not the Ends
                plt.plot(row[1],row[0], label='E', linewidth=1, color=(float(0), float(0), float(1)))
            else:
                plt.plot(row[1],row[0], label='E', linewidth=1, color=(float(c_r/100), float(c_g/100), float(0.5)))

                plt.plot(row[1][0], row[0][0], marker='.', markersize = "10" , color=(float(c_r/100), float(c_g/100), float(0.5)))
                plt.plot(row[1][-1], row[0][-1], marker='*', markersize = "10" , color=(float(c_r/100), float(c_g/100), float(0.5)))

        finally:
            s.close()

    plt.savefig("tra_reg_rgb" + str(RUN_ID) + "-" + str(random.randint(100, 999)) +".png")
    plt.show()



#mpl.use('macosx') #for the 3D magic

def stable_points_space_3d_rotate():

    fig = plt.figure(figsize=(300,300))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Total Abundance")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Abundance', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])

    index_A = 0

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        #print(len(row[1]))
        #print(len(row[0]))
        #print(len(rAx_prime[index_A]))


        row_abundance = []
        decompose = rAx_prime[index_A]
        #print(decompose)
        run_length = len(decompose[0])
        #print(run_length)


        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_abundance.append(current_sum)

        print(len(row_abundance))
        print(len(row[1]))
        print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_abundance[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], row_abundance[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_abundance[-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))

        index_A +=1

    #plt.savefig("3d_abundance_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

#stable_points_space_3d_rotate()


def stable_points_space_3d_rotate_number_alive():

    fig = plt.figure(figsize=(300,300))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Number Alive")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Alive', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])

    index_A = 0

    #print(rNumberAlive_prime)

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_alive = []

        decompose = rNumberAlive_prime[index_A]
        run_length = len(decompose[0])

        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_alive.append(current_sum)

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_alive[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_alive , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_alive, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], row_alive[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_alive[-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))

        index_A +=1

    #plt.savefig("3d_alive_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

#stable_points_space_3d_rotate_number_alive()




if __name__ == '__main__':
    #pool = Pool(processes=7)
    #pool.map(read_files_parallel, [_ for _ in data_archives])

    #def power_n(x, n):
    #    return x ** n

    #result = pool.starmap(power_n, [(x, 2) for x in range(20)])
    #print(result)

    #pool = Pool(processes=7)
    #stable_points_space_heat()

    #read_files_parallel(data_archives[0])
    stable_points_global_local()
    print("stable points global")
    #stable_points_global_local()
    print("abundance heat")
    #stable_points_space_final_abundance_rgb_heat()
    print("alive heat")
    #stable_points_space_final_alive_rgb_heat()

    ###### map_async call each function !


    #print("points space")
    #stable_points_space()
    #print("abundance heat map")
    #stable_points_space_final_abundance_rgb_heat()
    #print("alive heat map")
    #stable_points_space_final_alive_rgb_heat()

    heatmap = [[0 for _ in range(10)] for _ in range(10)]

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.show()
    for line in heatmap:
        print(line)

    # rE[0] = global
    # rE[1] = local
    # heatmap[int(row[1][0])][int(row[0][0])] = row_alive[-1]
    # lower
    # heatmap [global] [local]

    heatmap[1][5] = 10

    print("---")

    for line in heatmap:
        print(line)

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.show()