import sys
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from joblib import Parallel, delayed
import multiprocessing

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

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']

finally:
    s.close()

K = biotic_components_K
R = essential_range_R
P = external_perturbation_rate_P
start = time_start
end = time_end
step = 1
N = environment_components_N

ROUND = truncated_gaussian_ROUND
OEn = niche_width
OE = [OEn for _ in range(K)]
w = affects_w
u = optimum_condition_u
rNumberAlive = [[] for x in range(K)]
alpha = [[] for _ in range(K)]


alive_threshold = 0.5

#inputs = range(1000)
#def processInput(i):
#    return i * i
#num_cores = multiprocessing.cpu_count()
#results = Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)

def alive_species_count_no_gaussian():

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]


    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Alive Species Count - no Gaussian Truncation, AT = ' +  str(alive_threshold))
    plt.xlabel('EL')
    plt.ylabel('EG')
    #EG yAxis - EL xAxis
    #heatmap[20][30]=300
    #heatmap[50][70]=300
    Gindex = 0
    Lindex = 0


    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            ###################################################

            for each_species in range(biotic_components_K):
                abundance = 0
                if each_species in local_population_index:
                    abundance = (

                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2))))
                    )
                else:
                    abundance = (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))


                ###################################################

                if abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += 1

            Lindex += 1
        Lindex = 0
        Gindex += 1


    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', origin='lower')
    plt.savefig('alive_species_count_no_truncation_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()

    print("Alive non Truncated")

def alive_species_count_truncated_gaussian_red():

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    alive_threshold = 0

    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Alive Species Count - Truncated Gaussian (Red) with AT = ' +  str(alive_threshold))
    plt.xlabel('EL')
    plt.ylabel('EG')

    Gindex = 0
    Lindex = 0

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):
            print(Global, Local)
            ###################################################
            for each_species in range(biotic_components_K):
                abundance = 0

                if each_species in local_population_index:
                    abundance = (
                            round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                            *
                            round((math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                    )

                else:
                    abundance = round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                    #finding threshold with rounding ...


                            ###################################################

                if abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += 1

            Lindex += 1
        Lindex = 0
        Gindex += 1

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('alive_species_count_truncation_red_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()

    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Alive Threshold')
    plt.xlabel('EL')
    plt.ylabel('EG')

    for each_species in range(biotic_components_K):
        for _ in np.arange(0,R, 0.1):
            if(round((math.e) ** ((-1) * (((abs((_)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)>0):
                plt.axhline(y=(round((math.e) ** ((-1) * (((abs((_)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)), color='r', linestyle='-')
                break;

    plt.show()

def alive_species_count_truncated_gaussian_green():
    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]


    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Alive Species Count - * Truncation (Green) with AT = ' +  str(alive_threshold))
    plt.xlabel('EL')
    plt.ylabel('EG')

    Gindex = 0
    Lindex = 0

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            ###################################################
            for each_species in range(biotic_components_K):
                abundance = 0

                if each_species in local_population_index:
                    abundance = (
                        round(
                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            ,truncated_gaussian_ROUND)
                    )
                else:
                    abundance = round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)

                ###################################################

                if abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += 1

            Lindex += 1
        Lindex = 0
        Gindex += 1

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('alive_species_count_truncation_green_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()


def others():




    print("Alive Truncated")



    #round((math.e) ** ((-1) * (((abs((E[ai])-optimum_condition_u[ai][_])) ** 2) / (2*(OE[_]**2)))),truncated_gaussian_ROUND)

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Abundance of Species - no Gaussian Truncation with AT = ' +  str(alive_threshold))
    plt.xlabel('EL')
    plt.ylabel('EG')

    Gindex = 0
    Lindex = 0

    #print(optimum_condition_u)
    #print(OE)


    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            ###################################################
            sum_abundance = 0

            for each_species in range(biotic_components_K):
                abundance = 0
                if each_species in local_population_index:
                    abundance = (
                            round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                            *
                            round((math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                    )
                else:
                    abundance = round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)

                sum_abundance += abundance

                ###################################################

                if sum_abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += sum_abundance

            Lindex += 1
        Lindex = 0
        Gindex += 1

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', origin='lower')
    plt.savefig('abundance_of_species_no_truncation_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()

    print("Abundance non Truncated")


    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Abundance of Species - Truncated Gaussian (Red) with AT = '+  str(alive_threshold))
    plt.xlabel('EL')
    plt.ylabel('EG')

    Gindex = 0
    Lindex = 0

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            ###################################################
            sum_abundance = 0
            for each_species in range(biotic_components_K):
                abundance = 0

                if each_species in local_population_index:
                    abundance = (
                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2))))
                    )
                else:
                    abundance = (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                sum_abundance += abundance
                ###################################################

                #print(Gindex, Lindex, Global, Local, sum_abundance)
                if sum_abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += sum_abundance

            Lindex += 1
        Lindex = 0
        Gindex += 1

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('abundance_of_species_truncation_red_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()

    print("Abundance Truncated")



    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Abundance of Species - * Truncation (Green) with AT = ' +  str(alive_threshold))
    plt.xlabel('EL')
    plt.ylabel('EG')

    Gindex = 0
    Lindex = 0

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            ###################################################
            sum_abundance = 0
            for each_species in range(biotic_components_K):
                abundance = 0

                if each_species in local_population_index:
                    abundance = (
                        round(
                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            ,truncated_gaussian_ROUND)
                    )
                else:
                    abundance = (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                sum_abundance += abundance
                ###################################################

                #print(Gindex, Lindex, Global, Local, sum_abundance)
                if sum_abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += sum_abundance

            Lindex += 1
        Lindex = 0
        Gindex += 1

    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower')
    plt.savefig('abundance_of_species_truncation_green_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()



    Globals_L = []
    Locals_L = []
    Abundance_L = []
    Alives = []
    Alive_thres = []

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    Gindex = 0
    Lindex = 0

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            sum_abundance = 0
            num_alive = 0

            for each_species in range(biotic_components_K):
                abundance = 0

                if each_species in local_population_index:
                    abundance = (
                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2))))
                    )
                else:
                    abundance = (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                sum_abundance += abundance

                if sum_abundance > alive_threshold:
                    num_alive +=1

            Globals_L.append(Global)
            Locals_L.append(Local)
            Abundance_L.append(sum_abundance)
            Alives.append(num_alive)
            Alive_thres.append(alive_threshold)

            Lindex += 1
        Lindex = 0
        Gindex += 1


    fig = plt.figure(figsize=(20,20),dpi=300)
    ax = fig.add_subplot(321, projection='3d', adjustable='box')
    ax.set_title(label = "Abundance - no Truncation")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Abundance')
    ax.plot_trisurf(Locals_L, Globals_L, Alive_thres,cmap='inferno', edgecolor='none', alpha = 0.2);
    ax.plot_trisurf(Locals_L, Globals_L, Abundance_L,cmap='viridis', edgecolor='none');

    ax = fig.add_subplot(322, projection='3d', adjustable='box')
    ax.set_title(label = "Number Alive - no Truncation")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Alive')
    ax.plot_trisurf(Locals_L, Globals_L, Alive_thres,cmap='inferno', edgecolor='none', alpha = 0.2);
    ax.plot_trisurf(Locals_L, Globals_L, Alives,cmap='viridis', edgecolor='none');



    Globals_L = []
    Locals_L = []
    Abundance_L = []
    Alives = []
    Alive_thres = []

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    Gindex = 0
    Lindex = 0

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            sum_abundance = 0
            num_alive = 0

            ###################################################
            sum_abundance = 0

            for each_species in range(biotic_components_K):
                abundance = 0
                if each_species in local_population_index:
                    abundance = (
                            round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                            *
                            round((math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)
                    )
                else:
                    abundance = round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)

                sum_abundance += abundance

                ###################################################

                if sum_abundance > alive_threshold:
                    num_alive +=1

            Globals_L.append(Global)
            Locals_L.append(Local)
            Abundance_L.append(sum_abundance)
            Alives.append(num_alive)
            Alive_thres.append(alive_threshold)

            Lindex += 1
        Lindex = 0
        Gindex += 1


    ax = fig.add_subplot(323, projection='3d', adjustable='box')
    ax.set_title(label = "Abundance - Truncation (Red Example)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Abundance')
    ax.plot_trisurf(Locals_L, Globals_L, Alive_thres,cmap='inferno', edgecolor='none', alpha = 0.2);
    ax.plot_trisurf(Locals_L, Globals_L, Abundance_L,cmap='viridis', edgecolor='none');

    ax = fig.add_subplot(324, projection='3d', adjustable='box')
    ax.set_title(label = "Number Alive - Individually Truncated (Red Example)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Alive')
    ax.plot_trisurf(Locals_L, Globals_L, Alive_thres,cmap='inferno', edgecolor='none', alpha = 0.2);
    ax.plot_trisurf(Locals_L, Globals_L, Alives,cmap='viridis', edgecolor='none');


    Globals_L = []
    Locals_L = []
    Abundance_L = []
    Alives = []
    Alive_thres = []

    heatmap = [[0 for _ in np.arange(0,essential_range_R,step)] for _ in np.arange(0,essential_range_R,step)]

    Gindex = 0
    Lindex = 0

    for Global in np.arange(0,essential_range_R,step):
        for Local in np.arange(0,essential_range_R,step):

            sum_abundance = 0
            num_alive = 0

            ###################################################
            sum_abundance = 0

            for each_species in range(biotic_components_K):
                abundance = 0
                if each_species in local_population_index:
                    abundance = (
                        round(
                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(OE[each_species]**2))))
                            ,truncated_gaussian_ROUND)
                    )
                else:
                    abundance = round((math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(OE[each_species]**2)))),truncated_gaussian_ROUND)

                sum_abundance += abundance

                ###################################################

                if sum_abundance > alive_threshold:
                    num_alive +=1

            Globals_L.append(Global)
            Locals_L.append(Local)
            Abundance_L.append(sum_abundance)
            Alives.append(num_alive)
            Alive_thres.append(alive_threshold)

            Lindex += 1
        Lindex = 0
        Gindex += 1


    ax = fig.add_subplot(325, projection='3d', adjustable='box')
    ax.set_title(label = "Abundance - truncation (Green Example)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Abundance')
    ax.plot_trisurf(Locals_L, Globals_L, Alive_thres,cmap='inferno', edgecolor='none', alpha = 0.2);
    ax.plot_trisurf(Locals_L, Globals_L, Abundance_L,cmap='viridis', edgecolor='none');

    ax = fig.add_subplot(326, projection='3d', adjustable='box')
    ax.set_title(label = "Number Alive - truncation (Green Example)")
    ax.set_xlabel('X - EL')
    ax.set_ylabel('Y - EG')
    ax.set_zlabel('Z - Alive')
    ax.plot_trisurf(Locals_L, Globals_L, Alive_thres,cmap='inferno', edgecolor='none', alpha = 0.2);
    ax.plot_trisurf(Locals_L, Globals_L, Alives,cmap='viridis', edgecolor='none');




    print("Plotting ...")
    plt.savefig('z_alives_abundance_' +  str(truncated_gaussian_ROUND) + '_LPsz_'+ str(local_population_size) + '.png')
    plt.show()


if __name__ == '__main__':
    alive_species_count_no_gaussian()
    alive_species_count_truncated_gaussian_red()
    alive_species_count_truncated_gaussian_green()