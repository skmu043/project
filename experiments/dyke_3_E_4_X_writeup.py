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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
plt.rcParams["font.family"] = "Times New Roman"


# Generating ALL Parameters
SAMPLE_SIZE = 1
SAMPLE_STEP = 1
RUN_ID = int(time.time())

SPECIES_K   = 10000                  # ----------- Number of Biotic Components
RANGE_R     = 100                  # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 1                   # ----------- Time Step3
ENV_VARS    = 1                     # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0
ENV_START=[50]
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
global one
one = 0

Eg = ENV_START[0]

for s_i in range(SPECIES_K):
    a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

    if a_star < ALIVE_THRESHOLD:
        a_star = 0

    system_state[s_i] = a_star

for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START[_]

def rates_of_change_system_state(system_state):

    rate_of_change = system_state.copy()

    Eg = system_state[SPECIES_K+0]
    #ABUNDANCE
    for s_i in range(SPECIES_K):

        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

        rate_of_change[s_i] =  a_star - system_state[s_i]

    biotic_force_FG = 0
    #BIOTIC FORCE
    for s_i in range(SPECIES_K):
        biotic_force_FG += (system_state[s_i] * omega[0][s_i])
    rate_of_change[SPECIES_K+0] = (biotic_force_FG)

    if(add == 50):
        if(rate_of_change[SPECIES_K+0] < 30):
            rate_of_change[SPECIES_K+0] += 30
        if(rate_of_change[SPECIES_K+0] > 70):
            rate_of_change[SPECIES_K+0] -= 30
        if(rate_of_change[SPECIES_K+0] <= 70 or rate_of_change[SPECIES_K+0] >= 30 ):
            rate_of_change[SPECIES_K+0] -= 20


    return(rate_of_change)

if __name__ == '__main__':


    # system_state[0-99]  = species
    # system_state[100]   = environment
    # system_state[101]   = biotic_force
    # system_state[102]   = perturbing force
    # system_state[103]   = exponential temp runoff

    for y in range(5):

        E4 = []

        for x in range(4):
            one = x
            omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
            mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

            print(str(y)+"_"+str(x))
            #print(omega)
            #print(mu)
            f = open("omega-mu/"+str(y)+"_"+str(x)+".omega-mu_10", "w")
            f.write(str(omega))
            f.write("\n\n")
            f.write(str(mu))
            f.close()
            print("===")
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


            ALIVE_THRESHOLD=0
            results = [[] for _ in range(SPECIES_K+ENV_VARS)]
            times_steps=[]

            for step in np.arange(TIME_START, TIME_END, TIME_STEP):
                print(step)
                times_steps.append(step)
                for _ in range(SPECIES_K+ENV_VARS):
                    results[_].append(system_state[_])

                k1 = TIME_STEP * rates_of_change_system_state(system_state)
                k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
                k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
                add +=1
                k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)
                system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)
            ENV_VAR_ALIVE_ZERO_END = system_state[SPECIES_K+0]

            E4.append(results[-1])

        fig, ax = plt.subplots()

        plt.title(str(y)    , fontsize = 15)
        plt.xlabel('Time', fontsize = 15)
        plt.ylabel('Temperature', fontsize = 15)


        for entry in E4:
            ax.plot(times_steps,entry,'--',label = 'Perturbing Force')


        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        ax.tick_params(which='both', width=1)
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=4)

        plt.xlim([0, 100])
        #lt.ylim([-100, 50])

        #plt.legend(prop={'size': 12}, loc='lower right')
        plt.tight_layout()
        plt.savefig("omega-mu-jpg/"+str(y)+"_10.jpg")
        plt.show()


