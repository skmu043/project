import sys
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


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

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']
    temperatures = s['temperatures']
    time_step = s['time_step']

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

stable_points_space_heat()


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


def stable_points_space():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Regions', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)


    for row in rE_prime:

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

    #plt.savefig("tra_reg_rgb" + str(RUN_ID) + "-" + str(random.randint(100, 999)) +".png")
    plt.show()

stable_points_space()



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

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Final Abundance HeatMap', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    index_A = 0
    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_abundance = []
        decompose = rAx_prime[index_A]
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
        index_A +=1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')

    plt.show()

stable_points_space_final_abundance_rgb_heat()

def stable_points_space_final_alive_rgb_heat():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Final Alive HeatMap', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    index_A = 0
    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

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
        index_A +=1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='upper')

    plt.show()

stable_points_space_final_alive_rgb_heat()

mpl.use('macosx') #for the 3D magic

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

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

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

stable_points_space_3d_rotate()


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

stable_points_space_3d_rotate_number_alive()