import sys
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

if int(len(sys.argv)) != int(2):
    print("Args: shelve file name which contains all of > (K, R, P, E, start, end, step, EN, OE, LP_Z, RUN_ID)")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN=2, OE=5, LP_Z = (10 - 100), RUN_ID : epoch")
    print("exit")
    sys.exit()

s = shelve.open(str(sys.argv[1]))

try:
    exp_name            = s['exp_name']
    data_directory      = s['data_directory']
    shelve_file         = s['shelve_file']

    SPECIES_K           = s['SPECIES_K']
    RANGE_R             = s['RANGE_R']
    TIME_START          = s['TIME_START']
    TIME_END            = s['TIME_END']
    TIME_STEP           = s['TIME_STEP']
    ENV_VARS            = s['ENV_VARS']
    NICHE               = s['NICHE']
    SURVIVAL_THRESHOLD  = s['SURVIVAL_THRESHOLD']
    omega               = s['omega']
    mu                  = s['mu']
    ENV_START           = s['ENV_START']

finally:
    s.close()
# Initilize #

system_state                = np.zeros(SPECIES_K+ENV_VARS)
Eg                          = ENV_START

########################################################################################################################
# NICHE - Match Range on Niche Change

def fYaI(Xe, Ni, u, T):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if(abundance <= T):
        abundance = 0

    return(abundance)

def fXe(Ya, Ni, u):
    return (u + (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))))

def fXe_negative(Ya, Ni, u):
    return (u - (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))))

def fYaIx(Xe, Ni, u, NRange):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if((Xe >= u + NRange) or (Xe <= u - NRange)):
        abundance = 0

    return(abundance)

########################################################################################################################

if(NICHE == 5):
    for s_i in range(SPECIES_K):
        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
        if a_star < SURVIVAL_THRESHOLD:
            a_star = 0
        system_state[s_i] = a_star

###########=========================
if(NICHE == 10):
    for s_i in range(SPECIES_K):
        NRange = (fXe(SURVIVAL_THRESHOLD, 5, mu[0][s_i]) - mu[0][s_i])
        a_star = fYaIx(Eg, NICHE, mu[0][s_i], NRange)
        system_state[s_i] = a_star

#truncation_level = fYaI(fXe(0.2, 5, u), 20, u, 0))
###########=========================

########################################################################################################################
# Before env var gets added - the last spot is zero
abundance_start = sum(system_state)

number_alive_start = 0

for abundance_init in system_state:
    if abundance_init > 0:
        number_alive_start +=1


# Environment Init
for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START


def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change\
    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    rate_of_change = system_state.copy()

    Eg = system_state[SPECIES_K+0]

    if(NICHE == 5):
        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            if a_star < SURVIVAL_THRESHOLD:
                a_star = 0
            rate_of_change[s_i] =  a_star - system_state[s_i]
            #da/dt = a* - a
    if(NICHE == 10):
        for s_i in range(SPECIES_K):
            NRange = (fXe(SURVIVAL_THRESHOLD, 5, mu[0][s_i]) - mu[0][s_i])
            a_star = fYaIx(Eg, NICHE, mu[0][s_i], NRange)
            rate_of_change[s_i] =  a_star - system_state[s_i]
            #da/dt = a* - a

    biotic_force_FG = 0


    for s_i in range(SPECIES_K):
        # Global
        biotic_force_FG += (system_state[s_i] * omega[0][s_i])

    rate_of_change[SPECIES_K+0] = (biotic_force_FG)
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

    number_alive_end = 0

    aliveness = 0

    abundance_end = 0

    if(SURVIVAL_THRESHOLD > 0):
        aliveness = 0.000001

    for abundance_stream in results[:-1]:
        if abundance_stream[-1] > aliveness:
            number_alive_end +=1
            abundance_end += abundance_stream[-1]


    s = shelve.open(str(sys.argv[1]))


    try:
        s['NUMBER_ALIVE_START']     = number_alive_start
        s['NUMBER_ALIVE_END']       = number_alive_end
        s['ENV_END']                = results[SPECIES_K][-1]
        s['TOTAL_ABUNDANCE_START']  = abundance_start
        s['TOTAL_ABUNDANCE_END']    = abundance_end
    finally:
        s.close()

