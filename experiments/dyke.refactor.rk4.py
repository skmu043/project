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

    SPECIES_K = s['SPECIES_K']
    RANGE_R = s['RANGE_R']
    TIME_START = s['TIME_START']
    TIME_END = s['TIME_END']
    TIME_STEP = s['TIME_STEP']
    ENV_VARS = s['ENV_VARS']
    NICHE = s['NICHE']
    LOCAL_SIZE = s['LOCAL_SIZE']
    JI_ST = s['JI_ST']

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']

    omega = s['omega']
    mu = s['mu']
    local_population_index = s['local_population_index']

    ENV_START = s['ENV_START']

finally:
    s.close()

Eg = ENV_START[0]
El = ENV_START[1]

system_state = np.zeros(SPECIES_K+ENV_VARS)

#print(ENV_START)

#print(mu)
#print(omega)

################## INITIAL STATE
# Abundance Init


if(JI_ST==0): # JI Processing

    for s_i in range(SPECIES_K):
        if s_i in local_population_index:
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) \
                     * \
                     np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2))

            system_state[s_i] = a_star

        else :
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

            system_state[s_i] = a_star

if(JI_ST==0.2): # ST Processing

    for s_i in range(SPECIES_K):
        if s_i in local_population_index:
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) \
                     * \
                     np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2))

            if a_star <= 0.2:
                a_star = 0

            system_state[s_i] = a_star

        else :
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

            if a_star <= 0.2:
                a_star = 0

            system_state[s_i] = a_star

########################################################################################################################

# Environment Init
for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START[_]


def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change\
    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    rate_of_change = system_state.copy()

    Eg = system_state[SPECIES_K+0]
    El = system_state[SPECIES_K+1]

    if(JI_ST==0): # JI Processing
        for s_i in range(SPECIES_K):
            if s_i in local_population_index:
                a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) \
                         * \
                         np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2))

                rate_of_change[s_i] = a_star - system_state[s_i]

            else :
                a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

                rate_of_change[s_i] =  a_star - system_state[s_i]

    if(JI_ST==0.2): # ST Processing
        for s_i in range(SPECIES_K):
            if s_i in local_population_index:
                a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) \
                         * \
                         np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2))

                if a_star <= 0.2:
                    a_star = 0

                rate_of_change[s_i] = a_star - system_state[s_i]

            else :
                a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

                if a_star <= 0.2:
                    a_star = 0

                rate_of_change[s_i] =  a_star - system_state[s_i]


        #da/dt = a* - a
    biotic_force_FG = 0
    biotic_force_FL = 0

    for s_i in range(SPECIES_K):
        # Global
        biotic_force_FG += (system_state[s_i] * omega[0][s_i])
        # Local
        if s_i in local_population_index:
            biotic_force_FL += (system_state[s_i] * omega[1][s_i])

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


    s = shelve.open(shelve_file)


    try:
        #s['system_state'] = system_state
        s['results'] = results
        s['times_steps'] = times_steps
    finally:
            s.close()

