import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time

import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

#from numba import jit

# Generating ALL Parameters
SAMPLE_SIZE = 1
SAMPLE_STEP = 1
RUN_ID = int(time.time())

SPECIES_K   = 100                  # ----------- Number of Biotic Components
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


number_alive_global_start = 0
number_alive_start = 0

system_state = np.zeros(SPECIES_K+ENV_VARS)

Eg = ENV_START[0]

for s_i in range(SPECIES_K):

    a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

    if a_star < ALIVE_THRESHOLD:
        a_star = 0

    system_state[s_i] = a_star

    if a_star >= ALIVE_THRESHOLD:
        number_alive_global_start +=1


number_alive_start = number_alive_global_start

# Environment Init
for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START[_]

def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change\
    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    rate_of_change = system_state.copy()

    Eg = system_state[SPECIES_K+0]

    for s_i in range(SPECIES_K):

        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

        if a_star < ALIVE_THRESHOLD:
            a_star = 0

        rate_of_change[s_i] =  a_star - system_state[s_i]


        #da/dt = a* - a
    biotic_force_FG = 0

    for s_i in range(SPECIES_K):
        # Global
        biotic_force_FG += (system_state[s_i] * omega[0][s_i])

    rate_of_change[SPECIES_K+0] = (biotic_force_FG)

    #dE/dt = E* + F

    return(rate_of_change)



def plot_alphas():

    temperatures = []
    biotic_force = [[] for _ in range(SPECIES_K)]
    step = 0.01

    for x in np.arange (-50, RANGE_R+50, step):
        temperatures.append(x)

    for y in range(SPECIES_K):
        for x in np.arange (-50, RANGE_R+50, step):
            biotic_force[y].append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))) * omega[0][y])
            #biotic_force[y].append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force of 100 species with alive threshold : ' + str(ALIVE_THRESHOLD), fontsize=40)
    plt.xlabel('Environment Condition (temperature)', fontsize=40)
    plt.ylabel('Biotic Force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(SPECIES_K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)

    plt.show()

truncation = 0.2

def plot_alphas_truncated():


    temperatures = []
    alive_value = [[] for _ in range(SPECIES_K)]
    step = 0.01

    for x in np.arange (-50, RANGE_R+50, step):
        temperatures.append(x)

    for y in range(SPECIES_K):
        for x in np.arange (-50, RANGE_R+50, step):
            aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
            if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                alive_value[y].append(0)
            else:
                alive_value[y].append(aliveness * omega[0][y])


    plt.figure(figsize=(30,30))
    plt.title('Biotic Force of 100 species with alive threshold: ' + str(ALIVE_THRESHOLD), fontsize=40)
    plt.xlabel('Environment Condition (temperature)', fontsize=40)
    plt.ylabel('Biotic Force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(SPECIES_K):
        plt.plot(temperatures,alive_value[_])

    plt.plot(temperatures,np.sum((np.array(alive_value, dtype=float)), axis=0), lw=4)

    plt.show()



###################### ROOTS ########################################################
###################### ROOTS ########################################################


def f1(x):
    #return((x**3) + (2*(x**2)) - (2*x) - 5)
    #return(x**2 -1000)
    biotic_force = []
    for y in range(SPECIES_K):
        biotic_force.append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))) * omega[0][y])

    return(np.sum((np.array(biotic_force, dtype=float))))


    #print(xi," ",y[-1])

#TypeError: fsolve: there is a mismatch between the input and output shape of the 'func' argument 'f1'.Shape should be (2,) but it is (1,).

def plot_function():
    print("Plotting Sum  ... ")
    plt.figure(figsize=(20,10))
    plt.title('xy', fontsize=40)
    plt.xlabel('x', fontsize=40)
    plt.ylabel('y', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.axvline(x=0)
    plt.axhline(y=0)

    plt.plot(x,y, 'r-',label = 'roots')
    plt.show()


def plot_stable_points():

    x = []
    y = []

    X1 = -50
    Y1 = RANGE_R + 50

    for xi in np.arange(X1, Y1, 0.1):
        x.append(xi)
        y.append(f1(xi))

    print("Solving Roots ...")
    true_zeros = []

    for _ in range(RANGE_R):
        sol = optimize.root(f1, [_], jac=False, method='hybr')
        if(sol.x >=0 and sol.x <= RANGE_R):
            true_zeros.append(sol.x)


    #print("Solving ...")
    true_zeros = []
    sign_change = ""

    if(y[0] < 0):
        sign_change = "neg"
    if(y[0] > 0):
        sign_change = "pos"
    if(y[0] == 0):
        print("ZERO DETECTED")

    #print(sign_change)

    for _ in range(RANGE_R):
        sol = optimize.root(f1, [_], method='df-sane')
        if(sol.x >=0 and sol.x <= RANGE_R):
            true_zeros.append(sol.x)



    plt.figure(figsize=(20,10))
    plt.title('Roots', fontsize=40)
    plt.xlabel('temperature', fontsize=40)
    plt.ylabel('biotic force', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for stable in true_zeros:
        plt.axvline(x=stable)
    plt.axvline(x=0)
    plt.axhline(y=0)
    plt.plot(x,y, 'r-',label = 'biotic force')
    #plt.legend(loc=7, prop={'size': 30})
    plt.show()




#print(np.unique(np.array(true_zeros)))

###################### ROOTS ########################################################
###################### ROOTS ########################################################


###################### ROOTS ########################################################
###################### ROOTS ########################################################

# TRUNCATION


def f1_t(x):
    #return((x**3) + (2*(x**2)) - (2*x) - 5)
    #return(x**2 -1000)
    biotic_force = []
    for y in range(SPECIES_K):
        aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
        if(aliveness <= truncation and aliveness >= (-1 * truncation)):
            biotic_force.append(0 * omega[0][y])
        else:
            biotic_force.append(aliveness * omega[0][y])


        #biotic_force.append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))) * omega[0][y])

    return(np.sum((np.array(biotic_force, dtype=float))))


    #print(xi," ",y[-1])

#TypeError: fsolve: there is a mismatch between the input and output shape of the 'func' argument 'f1'.Shape should be (2,) but it is (1,).

def plot_function_t():
    print("Plotting Sum  ... ")
    plt.figure(figsize=(20,10))
    plt.title('xy', fontsize=40)
    plt.xlabel('x', fontsize=40)
    plt.ylabel('y', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.axvline(x=0)
    plt.axhline(y=0)

    plt.plot(x,y, 'r-',label = 'roots')
    plt.show()



def plot_stable_points_t():

    x = []
    y = []

    X1 = -50
    Y1 = RANGE_R + 50

    for xi in np.arange(X1, Y1, 0.1):
        x.append(xi)
        y.append(f1_t(xi))


    print("Solving Roots Truncated ...")
    true_zeros = []

    for _ in range(RANGE_R):
        sol = optimize.root(f1_t, [_], jac=False, method='hybr')
        if(sol.x >=0 and sol.x <= RANGE_R):
            true_zeros.append(sol.x)


    #print("Solving ...")
    true_zeros = []
    sign_change = ""

    if(y[0] < 0):
        sign_change = "neg"
    if(y[0] > 0):
        sign_change = "pos"
    if(y[0] == 0):
        print("ZERO DETECTED")

    #print(sign_change)

    for _ in range(RANGE_R):
        sol = optimize.root(f1_t, [_], method='df-sane')
        if(sol.x >=0 and sol.x <= RANGE_R):
            true_zeros.append(sol.x)

    plt.figure(figsize=(20,10))
    plt.title('Roots', fontsize=40)
    plt.xlabel('temperature', fontsize=40)
    plt.ylabel('biotic force', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for stable in true_zeros:
        plt.axvline(x=stable)
    plt.axvline(x=0)
    plt.axhline(y=0)
    plt.plot(x,y, 'r-',label = 'biotic force')
    #plt.legend(loc=7, prop={'size': 30})
    plt.show()




#print(np.unique(np.array(true_zeros)))

###################### ROOTS ########################################################


def plot_gaussian():

    ideal_temp = 50
    temp = []
    gaus = []
    for each_temp in np.arange(0,100,0.01):
        temp.append(each_temp)
        result = (math.e) ** ((-1) * (((abs(each_temp-ideal_temp)) ** 2) / (2*(NICHE**2))))
        gaus.append(result)

    plt.figure(figsize=(20,10))
    plt.title('The Gaussian Distribution', fontsize=40)
    plt.xlabel('Environment Condition (temperature)', fontsize=40)
    plt.ylabel('Alive Value', fontsize=40)
    plt.plot(temp,gaus, 'r-',label = 'The gaussian distribution')
    plt.show()

    print(gaus)

def gaus(each_temp):
    ideal_temp = 50
    temp = []
    gaus = []

    result = (math.e) ** ((-1) * (((abs(each_temp-ideal_temp)) ** 2) / (2*(NICHE**2))))

    if(result >= ALIVE_THRESHOLD):
        return(result)
    else:
        return(0)


def plot_gaussian_trunk():


    ideal_temp = 50

    temp = []
    gaus = []

    for each_temp in np.arange(0,100,0.01):
        temp.append(each_temp)
        result = (math.e) ** ((-1) * (((abs(each_temp-ideal_temp)) ** 2) / (2*(NICHE**2))))
        gaus.append(result)

    #plt.figure(figsize=(20,10))
    #plt.title('The Gaussian Distribution', fontsize=40)
    #plt.xlabel('Environment Condition (temperature)', fontsize=40)
    #plt.ylabel('Alive Value', fontsize=40)
    #plt.plot(temp,gaus, 'r-',label = 'The gaussian distribution')
    #plt.show()

    #plt.figure(figsize=(20,10))
    #plt.title('The Truncated Gaussian Distribution', fontsize=40)
    #plt.xlabel('Environment Condition (temperature)', fontsize=40)
    #plt.ylabel('Alive Value', fontsize=40)


    temp_t = []
    gaus_t = []
    for each_temp in np.arange(0,100,0.01):
        temp_t.append(each_temp)
        result = (math.e) ** ((-1) * (((abs(each_temp-ideal_temp)) ** 2) / (2*(NICHE**2))))
        if (result > ALIVE_THRESHOLD):
            gaus_t.append(result)
        else:
            gaus_t.append(0)

    #plt.axhline(y=ALIVE_THRESHOLD, color='g', linestyle='--')
    #plt.plot(temp_t,gaus,_t 'b-',label = 'The gaussian distribution')
    #plt.show()


    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('Abundance Values')
    #fig.set_size_inches(3, 1.5)
    ax1.plot(temp, gaus)
    ax1.set_title('Original')
    ax1.set_xlabel('Temperature')
    ax1.set_ylabel('Abundance')
    ax2.plot(temp_t, gaus_t)
    ax2.set_title('Survival Threshold')
    ax2.set_xlabel('Temperature')
    ax2.set_ylabel('Abundance')
    ax2.axhline(y=0.2, linestyle = '-')

    fig.show()





if __name__ == '__main__':

    ENV_VAR_ALIVE_ZERO_START = ENV_START[0]

    ENV_VAR_ALIVE_ONE_START = ENV_START[0]


    ALIVE_THRESHOLD=0

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


    ENV_VAR_ALIVE_ZERO_END = system_state[SPECIES_K+0]

    results_nt = results


    #================
    #plot_alphas()
    ALIVE_THRESHOLD=0.2
    #plot_alphas_truncated()
    #plot_stable_points()
    #plot_stable_points_t()
    #plot_gaussian()
    plot_gaussian_trunk()

    results = [[] for _ in range(SPECIES_K+ENV_VARS)]

    times_steps=[]

    system_state = np.zeros(SPECIES_K+ENV_VARS)

    Eg = ENV_START[0]

    for s_i in range(SPECIES_K):

        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

        if a_star < ALIVE_THRESHOLD:
            a_star = 0

        system_state[s_i] = a_star

        if a_star >= ALIVE_THRESHOLD:
            number_alive_global_start +=1

    for _ in range(ENV_VARS):
        system_state[SPECIES_K+_] = ENV_START[_]

    for step in np.arange(TIME_START, TIME_END, TIME_STEP):

        times_steps.append(step)

        for _ in range(SPECIES_K+ENV_VARS):
            results[_].append(system_state[_])

        k1 = TIME_STEP * rates_of_change_system_state(system_state)
        k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
        k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
        k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)

        system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)

    ENV_VAR_ALIVE_ONE_END = system_state[SPECIES_K+0]


    fig, axs = plt.subplots(2, 2, dpi=300, figsize=(10,20))
    fig.suptitle('Species Aliveness')
    myList = results_nt[:-1]
    for item in myList:
        axs[0,0].plot(times_steps,item)
    axs[0,0].set_title('Simulation with alive threshold')
    axs[0,0].set_xlabel('Time Steps')
    axs[0,0].set_ylabel('Alive Value for each species')

    myList = results[:-1]
    for item in myList:
        axs[0,1].plot(times_steps,item)

    axs[0,1].set_title('Simulation with alive threshold')
    axs[0,1].set_xlabel('Time Steps')
    axs[0,1].set_ylabel('Alive Value for each species')


    axs[1,0].set_title('Temperature')
    axs[1,0].set_xlabel('Time Steps')
    axs[1,0].set_ylabel('Temperature')
    axs[1,0].plot(times_steps,results_nt[-1])
    axs[1,0].plot(times_steps, results[-1])

    fig.show()




    number_alive_global_end = 0
    number_alive_end = 0

    for s_i in range(SPECIES_K):

        a_star = system_state[s_i]
        if a_star >= ALIVE_THRESHOLD:
            number_alive_global_end +=1

    number_alive_end = number_alive_global_end

