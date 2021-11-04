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
    local_population_index = s['local_population_index']

    global_start_temp = s['global_start_temp']
    local_start_temp = s['local_start_temp']

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
step = 0.001
N = environment_components_N
E = [global_start_temp, local_start_temp]
ROUND = truncated_gaussian_ROUND
OEn = niche_width
OE = [OEn for _ in range(K)]
w = affects_w
u = optimum_condition_u
rNumberAlive = [[] for x in range(K)]
alpha = [[] for _ in range(K)]


if __name__ == '__main__':



    #round((math.e) ** ((-1) * (((abs((E[ai])-optimum_condition_u[ai][_])) ** 2) / (2*(OE[_]**2)))),truncated_gaussian_ROUND)


    heatmap = [[0 for _ in np.arange(-50,essential_range_R+50,step)] for _ in np.arange(-50,essential_range_R+50,step)]

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



