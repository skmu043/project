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
SPECIES_K   = 100 # ----------- Number of Biotic Components
RANGE_R     = 100 # ----------- Essential Range
TIME_START  = 0   # ----------- Start of Simulation
TIME_END    = 200 # ----------- Length of Simulation
TIME_STEP   = 0.1 # ----------- Time Step
ENV_VARS    = 2   # ----------- Number of Environment Variables
NICHE_WIDTH = 5   # ----------- Niche Size
LOCAL_SIZE  = 0 # ----------- Local Population Size

En = []
for ei in range(ENV_VARS):
    En.append((random.uniform(10, RANGE_R)))
print(En)

affects_w = [[] for _ in range(ENV_VARS)]
for wi in range(ENV_VARS):
    affects_w[wi] = [random.uniform(-1, 1) for _ in range(SPECIES_K)]

optimum_condition_u = [[] for _ in range(ENV_VARS)]
for ui in range(ENV_VARS):
    optimum_condition_u[ui] = [random.uniform(0, RANGE_R) for _ in range(SPECIES_K)]

local_population_index = []
uniq_k = []
for x in range(int(LOCAL_SIZE/100 * SPECIES_K)):
    one = random.randint(0,SPECIES_K-1)
    while one in uniq_k:
        one = random.randint(0,SPECIES_K-1)
    uniq_k.append(one)
    local_population_index.append(one)
local_population_index.sort()

system_state = np.zeros(SPECIES_K+ENV_VARS)

plt.figure(figsize=(8,8), dpi=200)
plt.title('Gaussians')
plt.xlabel('Temperature')
plt.ylabel('Biotic Effect')

for E_index in range(ENV_VARS):
    for _ in range(SPECIES_K):
        tem = []
        abundance = []
        affects = []
        abundance_sum = 0

        for temp in np.arange(0,RANGE_R, 0.001):
            tem.append(temp)
            abundance.append(((math.e) ** ((-1) * (((abs((temp)-(optimum_condition_u[E_index][_]))) ** 2) / (2*(NICHE_WIDTH**2))))))
            affects.append(abundance[-1] * affects_w[0][_])
            abundance_sum += abundance[-1] * affects_w[0][_]
        plt.plot(tem,affects)

plt.show()

print(optimum_condition_u)
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

    #if step == 100:
    #    system_state[-1] += 10
    #    system_state[-2] -= 10


plt.figure(figsize=(8,8), dpi=200)
plt.title('Abundance')
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
