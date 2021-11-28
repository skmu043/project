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



system_state = np.zeros(SPECIES_K+ENV_VARS)



################## INITIAL STATE
# Abundance Init
for E_index in range(ENV_VARS):
    for _ in range(SPECIES_K):
        system_state[_] = ((math.e) ** ((-1) * (((abs((En[E_index])-(mu[E_index][_]))) ** 2) / (2*(NICHE_WIDTH**2)))))
# Environment Init
for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = En[_]
################## INITIAL STATE

print("System State : ", system_state)


def results_shelve():

    s = shelve.open(shelve_file)

    try:
        s['rAx'] = rAx
        s['rAxR'] = rAxR
        s['rNumberAlive'] = rNumberAlive
        s['alpha'] = alpha
        s['rF'] = rF
        s['rE'] = rE
        s['time'] = time
        s['OE'] = OE

    finally:
        s.close()


def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change\
    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    rate_of_change = system_state.copy()

    for s_i in range(SPECIES_K):
        if s_i in local_population_index:                     # Hard coded 0 and 1 -> see notes [E1], [E1, E2], [E3][E4][E5] - scale with Es like w, u there will be another
            Eg = system_state[SPECIES_K+0]
            El = system_state[SPECIES_K+1]

            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE_WIDTH ** 2 )) \
                     * \
                     np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE_WIDTH ** 2))

            if a_star < alive_threshold:
                a_star = 0

            rate_of_change[s_i] = a_star - system_state[s_i]

        else :
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE_WIDTH ** 2 ))

            if a_star < alive_threshold:
                a_star = 0

            rate_of_change[s_i] =  a_star - system_state[s_i]


        #da/dt = a* - a
    biotic_force_FG = 0
    biotic_force_FL = 0

    for _ in range(SPECIES_K):
        # Global
        biotic_force_FG += (system_state[_] * affects_w[0][_]) # >>> contains current rate of change for alpha
        # Local
        if _ in local_population_index:
            biotic_force_FL += (system_state[_] * affects_w[0][_]) # >>> contains current rate of change for alpha

    rate_of_change[SPECIES_K+0] = (biotic_force_FG)
    rate_of_change[SPECIES_K+1] = (biotic_force_FL)

    #dE/dt = E* + F

    return(rate_of_change)




if __name__ == '__main__':



    results = [[] for _ in range(SPECIES_K+ENV_VARS)]

    times_steps=[]

    for step in np.arange(TIME_START, TIME_END, TIME_STEP):

        times_steps.append(step)

        for _ in range(SPECIES_K+ENV_VARS):
            results[_].append(system_state[_])

        k1 = TIME_STEP * rates_of_change_system_state(system_state)
        k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
        k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
        k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)

        system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)









