import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
from matplotlib.gridspec import GridSpec

import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
#from numba import jit
plt.rcParams["font.family"] = "Times New Roman"


# Generating ALL Parameters
SAMPLE_SIZE = 1
SAMPLE_STEP = 1
RUN_ID = int(time.time())

SPECIES_K   = 100                  # ----------- Number of Biotic Components
RANGE_R     = 100                  # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 1                   # ----------- Time Step3
ENV_VARS    = 4                     # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0
ENV_START=[50,50,50,50]
omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

system_state = np.zeros(SPECIES_K+ENV_VARS)

# system_state[0-99]  = species
# system_state[100]   = environment
# system_state[101]   = biotic_force
# system_state[102]   = perturbing force
# system_state[103]   = exponential temp runoff
global biotic
biotic = []
global perturb
perturb = []
global exp_temp
exp_temp = []
global add
add = 0

Eg = ENV_START[0]

for s_i in range(SPECIES_K):
    a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

    system_state[s_i] = a_star

for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START[_]

def rates_of_change_system_state(system_state):

    rate_of_change = system_state.copy()

    E1 = system_state[SPECIES_K+0]
    E2 = system_state[SPECIES_K+0]
    E3 = system_state[SPECIES_K+0]
    E4 = system_state[SPECIES_K+0]

#ABUNDANCE
    for s_i in range(SPECIES_K):

        a_star1 = np.exp(- abs(E1-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
        a_star2 = np.exp(- abs(E2-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
        a_star3 = np.exp(- abs(E3-mu[2][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
        a_star4 = np.exp(- abs(E4-mu[3][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

        rate_of_change[s_i] =  a_star - system_state[s_i]

    biotic_force_FG_E1 = 0
    biotic_force_FG_E2 = 0
    biotic_force_FG_E3 = 0
    biotic_force_FG_E4 = 0
    #BIOTIC FORCE
    for s_i in range(SPECIES_K):
        biotic_force_FG_E1 += (system_state[s_i] * omega[0][s_i])
        biotic_force_FG_E2 += (system_state[s_i] * omega[1][s_i])
        biotic_force_FG_E3 += (system_state[s_i] * omega[2][s_i])
        biotic_force_FG_E4 += (system_state[s_i] * omega[3][s_i])

    rate_of_change[SPECIES_K+0] = (biotic_force_FG_E1)
    rate_of_change[SPECIES_K+1] = (biotic_force_FG_E2)
    rate_of_change[SPECIES_K+2] = (biotic_force_FG_E3)
    rate_of_change[SPECIES_K+3] = (biotic_force_FG_E4)

    if(add == 50):
        rate_of_change[SPECIES_K+0] = (biotic_force_FG_E1 + 10)
        rate_of_change[SPECIES_K+1] = (biotic_force_FG_E2 + 10)
        rate_of_change[SPECIES_K+2] = (biotic_force_FG_E3 + 10)
        rate_of_change[SPECIES_K+3] = (biotic_force_FG_E4 + 10)

    # if(add == 1):
    #     #BIOTIC
    #     biotic.append(biotic_force_FG * 10)
    #     #rate_of_change[SPECIES_K+0] = (biotic_force_FG + perturb[-1])
    #     #TEMP
    #     exp_temp.append(exp_temp[-1] + (exp_temp[-1] + perturb[-1])/20)
    #     #PERTURBING
    #     perturb.append(perturb[-1]+0.2)

    return(rate_of_change)

if __name__ == '__main__':


    # system_state[0-99]  = species
    # system_state[100]   = environment
    # system_state[101]   = biotic_force
    # system_state[102]   = perturbing force
    # system_state[103]   = exponential temp runoff

    for x in range(300):
        omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        # ======================================================================


        print(x)
        print(omega)
        print(mu)
        print("======================================================================")
        biotic = []
        perturb = []
        exp_temp = []
        add = 0



        system_state = np.zeros(SPECIES_K+ENV_VARS)

        Eg = ENV_START[0]

        # INIT ABUNDANCE
        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            system_state[s_i] = a_star
        #INIT ENV START
        for _ in range(ENV_VARS):
            system_state[SPECIES_K+_] = ENV_START[_]


        #INIT BIOTIC FORCE
        #biotic.append(0)
        #INIT PERTURBING FORCE
        perturb.append(0.2)
        #INIT EXPONENTIAL TEMP RUNOFF
        exp_temp.append(ENV_START[0])


        ALIVE_THRESHOLD=0
        results = [[] for _ in range(SPECIES_K+ENV_VARS)]
        times_steps=[]

        for step in np.arange(TIME_START, TIME_END, TIME_STEP):
            #print(step)
            times_steps.append(step)
            for _ in range(SPECIES_K+ENV_VARS):
                results[_].append(system_state[_])
            k1 = TIME_STEP * rates_of_change_system_state(system_state)
            k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
            k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
            k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)
            system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)
            add+=1
        ENV_VAR_ALIVE_ZERO_END = system_state[SPECIES_K+0]


        print(results)

        del perturb[-1]
        #del biotic[-1]
        del exp_temp[-1]




        plt.figure(figsize=(20,10))
        plt.title(str(x), fontsize=40)
        plt.xlabel('Time Steps', fontsize=40)
        plt.ylabel('Temperature, Perturb', fontsize=40)
        # index = 0
        # for item in results:
        #     plot_it = []
        #     if(index <=99) :
        #         for thing in item:
        #             plot_it.append(thing * 100)
        #         #plt.plot(times_steps,plot_it)
        #     else:
        #         plt.plot(times_steps,results[-1], 'r-',label = 'System Temperature', linewidth=5)
        #         plt.plot(times_steps,results[-2], 'g-',label = 'System Temperature', linewidth=5)
        #         plt.plot(times_steps,results[-3], 'b-',label = 'System Temperature', linewidth=5)
        #         plt.plot(times_steps,results[-4], 'k-',label = 'System Temperature', linewidth=5)
        #     index+=1
        #
        plt.plot(times_steps,results[-1], 'r-',label = 'System Temperature', linewidth=5)
        plt.plot(times_steps,results[-2], 'g-',label = 'System Temperature', linewidth=5)
        plt.plot(times_steps,results[-3], 'b-',label = 'System Temperature', linewidth=5)
        plt.plot(times_steps,results[-4], 'k-',label = 'System Temperature', linewidth=5)
        #plt.xlim([0, 200])
        #plt.ylim([-100, 150])
        plt.legend()
        plt.tight_layout()
        plt.savefig(str(x)+".jpg")
        #plt.show()

