import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import time

x = np.linspace(1,2,2)
y = np.linspace(3,4,2)
xx,yy = np.meshgrid(x,y)
result = xx * 2 + yy * 5

SAMPLE_SIZE = 1
SAMPLE_STEP = 1
RUN_ID      = time.time()
SPECIES_K   = 10                   # ----------- Number of Biotic Components
RANGE_R     = 100                   # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 0.1                   # ----------- Time Step
ENV_VARS    = 2                     # ----------- Number of Environment Variables
NICHE_WIDTH = 5                     # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)

En = [(random.uniform(10, RANGE_R)) for env_start_temp in range(ENV_VARS)]

affects_w = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

optimum_condition_u = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

local_population_index = []
uniq_k = []
for _ in range(int(LOCAL_SIZE/100 * SPECIES_K)):
    local_species = random.randint(0,SPECIES_K-1)
    while local_species in local_population_index:
        local_species = random.randint(0,SPECIES_K-1)
    local_population_index.append(local_species)
local_population_index.sort()

system_state = np.zeros(SPECIES_K+ENV_VARS)

plt.figure(figsize=(8,8), dpi=200)
plt.title('Biotic Effect of Each Species (a*w)')
plt.xlabel('Temperature')
plt.ylabel('Biotic Effect')

for E_index in range(ENV_VARS):
    for _ in range(SPECIES_K):
        tem = []
        abundance = []
        affects = []

        for temp in np.arange(0,RANGE_R, 0.001):
            tem.append(temp)
            abundance.append(((math.e) ** ((-1) * (((abs((temp)-(optimum_condition_u[E_index][_]))) ** 2) / (2*(NICHE_WIDTH**2))))))
            affects.append(abundance[-1] * affects_w[0][_])
        plt.plot(tem,affects)

plt.show()

################## INITIAL STATE
# Abundance Init
for E_index in range(ENV_VARS):
    for _ in range(SPECIES_K):
        system_state[_] = ((math.e) ** ((-1) * (((abs((En[E_index])-(optimum_condition_u[E_index][_]))) ** 2) / (2*(NICHE_WIDTH**2)))))
# Environment Init
for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = En[_]
################## INITIAL STATE

print("System State : ", system_state)

def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change\
    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    new_system_state = system_state.copy()

    for _ in range(SPECIES_K):
        if _ in local_population_index:                     # Hard coded 0 and 1 -> see notes [E1], [E1, E2], [E3][E4][E5] - scale with Es like w, u there will be another
            new_system_state[_] =  ((math.e) ** ((-1) * (((abs((system_state[SPECIES_K+0])-(optimum_condition_u[0][_]))) ** 2) / (2*(NICHE_WIDTH**2))))) \
                                   *\
                                   ((math.e) ** ((-1) * (((abs((system_state[SPECIES_K+1])-(optimum_condition_u[1][_]))) ** 2) / (2*(NICHE_WIDTH**2)))))\
                                   - system_state[_]
        else :
            new_system_state[_] =  ((math.e) ** ((-1) * (((abs((system_state[SPECIES_K+0])-(optimum_condition_u[0][_]))) ** 2) / (2*(NICHE_WIDTH**2))))) - system_state[_]

        #da/dt = a* - a
    biotic_force_FG = 0
    biotic_force_FL = 0

    for _ in range(SPECIES_K):
        # Global
        biotic_force_FG += (system_state[_] * affects_w[0][_]) # >>> contains current rate of change for alpha
        # Local
        if _ in local_population_index:
            biotic_force_FL += (system_state[_] * affects_w[0][_]) # >>> contains current rate of change for alpha

    new_system_state[SPECIES_K+0] = (biotic_force_FG)
    new_system_state[SPECIES_K+1] = (biotic_force_FL)

    #dE/dt = E* + F

    return(new_system_state)

results = [[] for _ in range(SPECIES_K+ENV_VARS)]

times_steps=[]

for step in np.arange(TIME_START, TIME_END, TIME_STEP):

    times_steps.append(step)

    for _ in range(SPECIES_K+ENV_VARS):
        results[_].append(system_state[_])

    k1 = TIME_STEP * (rates_of_change_system_state(system_state))
    k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
    k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
    k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)

    system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)


plt.figure(figsize=(8,8), dpi=200)
plt.title('Abundance for each Species')
plt.xlabel('Time Steps')
plt.ylabel('abundance values')

for _ in range(SPECIES_K):
    if _ in local_population_index:
        plt.plot(times_steps, results[_], 'r-')
    else:
        plt.plot(times_steps, results[_], 'b-')

plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.title('Environment Variables')
plt.xlabel('Time Steps')
plt.ylabel('Temp/PH')
plt.plot(times_steps, results[-2], 'k-', label = 'global')
plt.plot(times_steps, results[-1], 'b-', label = 'local')
plt.legend()
plt.show()


alive_threshold = 0.1

def alive_species_count_no_at_truncation_start():

    heatmap = [[0 for _ in np.arange(0,RANGE_R,TIME_STEP)] for _ in np.arange(0,RANGE_R,TIME_STEP)]

    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Alive Species Count : with AT Truncation, AT = ' +  str(alive_threshold))
    plt.xlabel('EL - PH')
    plt.ylabel('EG - Temp')

    xticks = np.round(np.linspace(0,RANGE_R/TIME_STEP, 7), 1)
    xtick_labels = [str(round(i)) for i in np.linspace(0, RANGE_R, 7)]
    plt.xticks(xticks, xtick_labels)

    yticks = np.round(np.linspace(0,RANGE_R/TIME_STEP, 7), 1)
    ytick_labels = [str(round(i)) for i in np.linspace(0, RANGE_R, 7)]
    plt.yticks(yticks, ytick_labels)

    Gindex = 0
    Lindex = 0


    for Global in np.arange(0,RANGE_R,TIME_STEP):
        for Local in np.arange(0,RANGE_R,TIME_STEP):

            ###################################################
            for each_species in range(SPECIES_K):
                abundance = 0
                if each_species in local_population_index:
                    abundance = (

                            (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(NICHE_WIDTH**2))))
                            *
                            (math.e) ** ((-1) * (((abs((Local)-optimum_condition_u[1][each_species])) ** 2) / (2*(NICHE_WIDTH**2))))
                    )
                else:
                    abundance = (math.e) ** ((-1) * (((abs((Global)-optimum_condition_u[0][each_species])) ** 2) / (2*(NICHE_WIDTH**2))))
            ###################################################

                if abundance > alive_threshold:
                    heatmap[Gindex][Lindex] += 1

            Lindex += 1
        Lindex = 0
        Gindex += 1


    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', origin='lower')
    plt.show()

    print("Alive non Truncated")

if __name__ == '__main__':
    alive_species_count_no_at_truncation_start()
