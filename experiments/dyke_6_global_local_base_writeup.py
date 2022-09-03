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
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

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
TIME_STEP   = 1                   # ----------- Time Step3                  # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 20                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0.00001
SURVIVAL_THRESHOLD = 0 # First run is 0 >> JI , second run updates before run to 0.2 for ST
ENV_START=[50,50] # Global E then Local E
ENV_VARS    = len(ENV_START)
omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

number_alive_global_start = 0
number_alive_start = 0

system_state = np.zeros(SPECIES_K+ENV_VARS)

Eg = ENV_START[0]
El = ENV_START[1]

# Select Local Sub Population

local_population_size = int(LOCAL_SIZE/100 * SPECIES_K)
local_population_index = []
local_population_index = random.sample(range(SPECIES_K), local_population_size)

for s_i in range(SPECIES_K):

    a_star = 0
    a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

    if(SURVIVAL_THRESHOLD == 0.2):
        if(a_star <= 0.2):
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
    El - system_state[SPECIES_K+1]

    for s_i in range(SPECIES_K):

        a_star = 0
        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * np.exp(- abs(El-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

        if(SURVIVAL_THRESHOLD == 0.2):
            if(a_star <= 0.2):
                a_star = 0

        system_state[s_i] = a_star

        rate_of_change[s_i] =  a_star - system_state[s_i]

    biotic_force_E1 = 0
    biotic_force_E2 = 0


    #each species has two affects, one for each ENV variable

    for s_i in range(SPECIES_K):
        biotic_force_E1 += (system_state[s_i] * omega[0][s_i]) # Effects on E1
    for s_i in range(SPECIES_K):
        biotic_force_E2 += (system_state[s_i] * omega[1][s_i]) # Effects on E2

    rate_of_change[SPECIES_K+0] = (biotic_force_E1)
    rate_of_change[SPECIES_K+1] = (biotic_force_E2)
    #dE/dt = E* + F

    return(rate_of_change)


def st_trunc(x):
    if x<=0.2:
        return 0
    return x

def plot_JI_ST_Biotic_ALL_E():

    E1 = np.arange(0,100,0.1)
    E2 = np.arange(0,100,0.1)

    E1, E2 = np.meshgrid(E1, E2)

    biotic_force_global = (math.e) ** ((-1) * (((abs(E1-500)) ** 2) / (2*(NICHE**2)))) + (E2*0) *0

    for s_i in range(SPECIES_K):
        E1_biotic = np.exp(- abs(E1-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * omega[0][s_i] + (E2*0)
        E2_biotic = np.exp(- abs(E2-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * omega[1][s_i] + (E1*0)

        biotic_force_global += E1_biotic
        biotic_force_global += E2_biotic


    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('JI Model',fontsize = 14)
    ax.set_xlabel('E1')
    ax.set_ylabel('E2')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(E1, E2, biotic_force_global, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()




    E1 = np.arange(0,100,0.1)
    E2 = np.arange(0,100,0.1)

    E1, E2 = np.meshgrid(E1, E2)

    biotic_force_global = (math.e) ** ((-1) * (((abs(E1-500)) ** 2) / (2*(NICHE**2)))) + (E2*0) *0

    for s_i in range(SPECIES_K):
        E1_biotic = np.exp(- abs(E1-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) + (E2*0)
        E2_biotic = np.exp(- abs(E2-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 )) + (E1*0)

        #print(E1_biotic.shape)
        index1=0
        index2=0
        for item in E1_biotic:
            for thing in item:
                if thing <=0.2:
                    E1_biotic[index1][index2] = 0
                index2+=1
            index1+=1
            index2=0

       # print(index1, index2)


        E1_biotic = E1_biotic * omega[0][s_i]
        E2_biotic = E2_biotic * omega[1][s_i]

        biotic_force_global += E1_biotic
        biotic_force_global += E2_biotic


    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('ST Model',fontsize = 14)
    ax.set_xlabel('E1')
    ax.set_ylabel('E2')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(E1, E2, biotic_force_global, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

plot_JI_ST_Biotic_ALL_E()


def plot_attractor_global_local_mesh():

    ideal_temp_global = 50
    ideal_temp_local = 80

    global_temp = np.arange(0,100,0.01)
    local_temp = np.arange(0,100,0.01)

    global_temp, local_temp = np.meshgrid(global_temp, local_temp)

    abundance_global = (math.e) ** ((-1) * (((abs(global_temp-500)) ** 2) / (2*(NICHE**2)))) + (local_temp*0)
    abundance_local = (math.e) ** ((-1) * (((abs(global_temp-500)) ** 2) / (2*(NICHE**2)))) * (math.e) ** ((-1) * (((abs(local_temp-500)) ** 2) / (2*(NICHE**2))))

    for s_i in range(SPECIES_K):
        if s_i in local_population_index:
            abundance = np.exp(- abs(global_temp-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * np.exp(- abs(local_temp-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            abundance_local += abundance
        if s_i not in local_population_index:
            abundance = np.exp(- abs(global_temp-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) + (local_temp*0)
            abundance_global += abundance


    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('Global Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(global_temp, local_temp, abundance_global, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('Local Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(global_temp, local_temp, abundance_local, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('Combined Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(global_temp, local_temp, abundance_global + abundance_local, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()



def plot_attractor_global_local_mesh_effects():


    global_temp = np.arange(0,100,0.01)
    local_temp = np.arange(0,100,0.01)

    global_temp, local_temp = np.meshgrid(global_temp, local_temp)

    biotic_force_global = (math.e) ** ((-1) * (((abs(global_temp-500)) ** 2) / (2*(NICHE**2)))) + (local_temp*0)
    biotic_force_local = (math.e) ** ((-1) * (((abs(global_temp-500)) ** 2) / (2*(NICHE**2)))) * (math.e) ** ((-1) * (((abs(local_temp-500)) ** 2) / (2*(NICHE**2))))

    for s_i in range(SPECIES_K):
        if s_i in local_population_index:
            biotic_force = np.exp(- abs(global_temp-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * np.exp(- abs(local_temp-mu[1][s_i]) ** 2 / ( 2 * NICHE ** 2 )) * omega[1][s_i]
            biotic_force_local += biotic_force
        if s_i not in local_population_index:
            biotic_force = np.exp(- abs(global_temp-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 )) + (local_temp*0) * omega[0][s_i]
            biotic_force_global += biotic_force


    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('Global Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(global_temp, local_temp, biotic_force_global, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('Local Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(global_temp, local_temp, biotic_force_local, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

    fig = plt.figure(dpi=300, figsize=(10,10))
    ax = fig.gca(projection='3d')
    ax.set_title('Combined Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    #plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    surf = ax.plot_surface(global_temp, local_temp, biotic_force_global + biotic_force_local, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    #ax.set_zlim(-1.01, 1.01)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

plot_attractor_global_local_mesh_effects()
plot_attractor_global_local_mesh()

def plot_abundance_global_at_50_optimal():

    ideal_temp = 50
    temp = []
    gaus = []
    for each_temp in np.arange(0,100,0.01):
        temp.append(each_temp)
        result = (math.e) ** ((-1) * (((abs(each_temp-ideal_temp)) ** 2) / (2*(NICHE**2))))
        gaus.append(result)

    plt.figure(figsize=(20,10))
    plt.title('Global Population', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('Abundance', fontsize=40)
    plt.plot(temp,gaus, 'b',label = 'The gaussian distribution')
    plt.show()


#plot_abundance_global_at_50_optimal()

def plot_abundance_local_at_50_optimal_mesh():

    ideal_temp_global = 50
    ideal_temp_local = 80
    fig, ax = plt.subplots(1, dpi=300, figsize=(20,20))
    global_temp = np.arange(0,100,0.01)
    local_temp = np.arange(0,100,0.01)
    global_temp, local_temp = np.meshgrid(global_temp, local_temp)
    abundance = (math.e) ** ((-1) * (((abs(global_temp-ideal_temp_global)) ** 2) / (2*(NICHE**2)))) * (math.e) ** ((-1) * (((abs(local_temp-ideal_temp_local)) ** 2) / (2*(NICHE**2))))
    ax.set_title('Local Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    plt.pcolormesh(global_temp, local_temp, abundance, cmap='viridis')
    plt.colorbar()
    plt.show()


plot_abundance_local_at_50_optimal_mesh()

def plot_abundance_local_at_50_optimal():

    ideal_temp_global = 50
    ideal_temp_local = 80
    fig, ax = plt.subplots(1, dpi=300, figsize=(20,20), subplot_kw={"projection": "3d"})
    global_temp = np.arange(0,100,0.01)
    local_temp = np.arange(0,100,0.01)
    global_temp, local_temp = np.meshgrid(global_temp, local_temp)
    abundance = (math.e) ** ((-1) * (((abs(global_temp-ideal_temp_global)) ** 2) / (2*(NICHE**2)))) * (math.e) ** ((-1) * (((abs(local_temp-ideal_temp_local)) ** 2) / (2*(NICHE**2))))
    ax.set_title('Local Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')

    surf = ax.plot_surface(global_temp, local_temp, abundance, cmap=plt.get_cmap('plasma'), linewidth=0, antialiased=False)
    ax.zaxis.set_major_formatter('{x:.02f}')
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()

plot_abundance_local_at_50_optimal()

def plot_20_local_population():

    fig, ax = plt.subplots(1, dpi=300, figsize=(20,20))
    global_temp = np.arange(0,100,0.01)
    local_temp = np.arange(0,100,0.01)
    global_temp, local_temp = np.meshgrid(global_temp, local_temp)
    abundance = (math.e) ** ((-1) * (((abs(global_temp-mu[0][local_population_index[0]])) ** 2) / (2*(NICHE**2)))) * (math.e) ** ((-1) * (((abs(local_temp-mu[1][local_population_index[0]])) ** 2) / (2*(NICHE**2))))

    local_population_index_first_removed = local_population_index.copy()
    del local_population_index_first_removed[0]

    abundances = []
    abundances.append(abundance)

    for s_i in range(SPECIES_K):
        if s_i in local_population_index_first_removed:
            abundance = (math.e) ** ((-1) * (((abs(global_temp-mu[0][s_i])) ** 2) / (2*(NICHE**2)))) * (math.e) ** ((-1) * (((abs(local_temp-mu[1][s_i])) ** 2) / (2*(NICHE**2))))
            abundances.append(abundance)


    ax.set_title('Local Population',fontsize = 14)
    ax.set_xlabel('Global Temperature')
    ax.set_ylabel('Local Temperature')
    for each_abundance in abundances:
        plt.pcolormesh(global_temp, local_temp, each_abundance, cmap='viridis')
    plt.colorbar()
    plt.show()

plot_20_local_population()

if __name__ == '__main__':


    omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
    mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

    ji_alives_at_start = [0 for _ in range(SPECIES_K)]
    ji_alives_at_end   = [0 for _ in range(SPECIES_K)]
    st_alives_at_start = [0 for _ in range(SPECIES_K)]
    st_alives_at_end   = [0 for _ in range(SPECIES_K)]

    system_state = np.zeros(SPECIES_K+ENV_VARS)

    Eg = ENV_START[0]

    for _ in range(ENV_VARS):
        system_state[SPECIES_K+_] = ENV_START[_]

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
    #plot_gaussian_trunk()
    #plot_temps()
    SURVIVAL_THRESHOLD=0.2
    #plot_alphas_truncated()

    results = [[] for _ in range(SPECIES_K+ENV_VARS)]
    times_steps=[]
    system_state = np.zeros(SPECIES_K+ENV_VARS)
    Eg = ENV_START[0]

    for s_i in range(SPECIES_K):
        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
        if a_star < ALIVE_THRESHOLD:
            a_star = 0
        system_state[s_i] = a_star

        if(a_star > 0.00001):
            st_alives_at_start[s_i] = 1

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

    print(results_nt[-1][-1])
    print(results[-1][-1])

    ji_final_results = results_nt

    for s_i in range(SPECIES_K):
        if(ji_final_results[s_i][-1] > 0):
            ji_alives_at_end[s_i] = 1

    print(results[-1])
    st_final_results = results
    for s_i in range(SPECIES_K):
        if(st_final_results[s_i][-1] > 0.00001):
            st_alives_at_end[s_i] = 1

    print(ji_alives_at_start)
    print(ji_alives_at_end)
    print(st_alives_at_start)
    print(st_alives_at_end)

    temperatures = []
    biotic_force = [[] for _ in range(SPECIES_K)]
    ji_alive_start_till_end = []
    ji_alive_only_start = []
    ji_not_alive_start_alive_end = []

    st_alive_start_till_end = []
    st_alive_only_start = []
    st_not_alive_start_alive_end = []

    step = 0.01

    for x in np.arange (-25, RANGE_R+25, step):
        temperatures.append(x)

    for y in range(SPECIES_K):
        if(ji_alives_at_start[y] == 1 and ji_alives_at_end[y]==1): # alive start and made it to end
            data = []
            for x in np.arange (-25, RANGE_R+25, step):
                data.append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
            ji_alive_start_till_end.append(data)
    for y in range(SPECIES_K):
        if(ji_alives_at_start[y] == 1 and ji_alives_at_end[y]==0): # alive start did not make it to end
            data = []
            for x in np.arange (-25, RANGE_R+25, step):
                data.append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
            ji_alive_only_start.append(data)
    for y in range(SPECIES_K):
        if(ji_alives_at_start[y] == 0 and ji_alives_at_end[y]==1): # not alive but alive at end
            data = []
            for x in np.arange (-25, RANGE_R+25, step):
                data.append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
            ji_not_alive_start_alive_end.append(data)




    for y in range(SPECIES_K):
        if(st_alives_at_start[y] == 1 and st_alives_at_end[y]==1): # alive start and made it to end
            data = []
            for x in np.arange (-25, RANGE_R+25, step):
                aliveness = ((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
                if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                    data.append(0)
                else:
                    data.append(aliveness)
            st_alive_start_till_end.append(data)
    for y in range(SPECIES_K):
        if(st_alives_at_start[y] == 1 and st_alives_at_end[y]==0): # alive start did not make it to end
            data = []
            for x in np.arange (-25, RANGE_R+25, step):
                aliveness = ((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
                if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                    data.append(0)
                else:
                    data.append(aliveness)
            st_alive_only_start.append(data)
    for y in range(SPECIES_K):
        if(st_alives_at_start[y] == 0 and st_alives_at_end[y]==1): # not alive but alive at end
            data = []
            for x in np.arange (-25, RANGE_R+25, step):
                aliveness = ((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
                if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                    data.append(0)
                else:
                    data.append(aliveness)
            st_not_alive_start_alive_end.append(data)



    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    #fig.suptitle('Abundance for 20 species', fontsize=30)
    #fig.set_size_inches(3, 1.5)

    ax1.set_title('JI Model', fontsize=35)
    ax1.set_xlabel('Temperature', fontsize=30)
    ax1.set_ylabel('Abundance', fontsize=30)
    for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
        label.set_fontsize(23)

    for item in ji_alive_start_till_end:
        ax1.plot(temperatures,item, 'b', label='Alive from start till end')
    for item in ji_alive_only_start:
        ax1.plot(temperatures,item, 'r', label='Alive at start but not alive at end')
    for item in ji_not_alive_start_alive_end:
        ax1.plot(temperatures,item, 'g', label='Not alive at start but alive at end')

    ax2.set_title('ST Model', fontsize=35)
    ax2.set_xlabel('Temperature', fontsize=30)
    ax2.set_ylabel('Abundance', fontsize=30)
    for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
        label.set_fontsize(23)

    for item in st_alive_start_till_end:
        ax2.plot(temperatures,item, 'b',label='Alive from start till end')
    for item in st_alive_only_start:
        ax2.plot(temperatures,item, 'r', label='Alive at start but not alive at end')
    for item in st_not_alive_start_alive_end:
        ax2.plot(temperatures,item, 'g', label='Not alive at start but alive at end')
    fig.tight_layout()
    #fig.legend()
    fig.show()


#if((results_nt[-1][-1] > 110 or results_nt[-1][-1] < -5) and (results[-1][-1] > 110 or results[-1][-1] < -5)):
    if(1):
        #or
       # ((results_nt[-1][-1] < 100 and results_nt[-1][-1] > 0) and (results[-1][-1] > 100 or results[-1][-1] < 0))):

    #print(omega)
    #print(mu)

    # GHOST NUMBERS

    #if((results[-1][-1] > 106 or results[-1][-1] < -6)):

        print(omega)
        print(mu)
        print(Eg)

        for line in results:
            print(line[-2])

        print("=================================================")

        fig = plt.figure(dpi=300, figsize=(20,10))
        #fig.suptitle('Species Aliveness ' + str(sim))
        #fig.suptitle('A simulation run with 100 biotic components', fontsize=20)

        gs = fig.add_gridspec(2,2)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[1, :])

        SIZE = 17



        myList = results_nt[:-1]
        for item in myList:
            ax1.plot(times_steps,item)
        ax1.set_title('JI Model', fontsize=20)
        ax1.set_xlabel('Time Steps', fontsize=19)
        ax1.set_ylabel('Abundance', fontsize=19)
        for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
            label.set_fontsize(SIZE)
        #ax1.set_ylim([0, 1])
        myList = results[:-1]
        for item in myList:
            ax2.plot(times_steps,item)
        ax2.set_title('ST Model', fontsize=20)
        ax2.set_xlabel('Time Steps', fontsize=19)
        ax2.set_ylabel('Abundance', fontsize=19)
        for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
            label.set_fontsize(SIZE)
        #ax2.set_ylim([0, 1])
        ax3.set_title('The Environment Condition',fontsize=20)
        ax3.set_xlabel('Time Steps', fontsize=19)
        ax3.set_ylabel('Temperature', fontsize=19)
        ax3.plot(times_steps,results_nt[-1], "b", label = "JI Model")
        ax3.plot(times_steps, results[-1],"k", label = "ST Model")
        for label in (ax3.get_xticklabels() + ax3.get_yticklabels()):
            label.set_fontsize(SIZE)

        #ax3.set_ylim([0, 100])
        plt.subplots_adjust(hspace=0.1)
        ax3.legend(prop={'size': 15})
        fig.tight_layout()
        fig.show()


    #number_alive_global_end = 0
        #number_alive_end = 0

        #for s_i in range(SPECIES_K):

        #    a_star = system_state[s_i]
        #    if a_star >= ALIVE_THRESHOLD:
        #        number_alive_global_end +=1

        #number_alive_end = number_alive_global_end

