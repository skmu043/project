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

SAMPLE_SIZE = int(1)
SAMPLE_STEP = int(1)
RUN_ID = int(time.time())
biotic_components_K = int(100)
essential_range_R = int(100)
external_perturbation_rate_P = int(0)
time_start = int(0)
time_end = int(200)
time_step = float(1)
environment_components_N = int(2)
niche_width = int(5)
local_population_size = int(40)

En = []
for ei in range(environment_components_N):
    En.append((random.uniform(10, essential_range_R)))

affects_w = [[] for _ in range(environment_components_N)]
for wi in range(environment_components_N):
    affects_w[wi] = [random.uniform(-1, 1) for _ in range(biotic_components_K)]

optimum_condition_u = [[] for _ in range(environment_components_N)]
for ui in range(environment_components_N):
    optimum_condition_u[ui] = [random.uniform(0, essential_range_R) for _ in range(biotic_components_K)]

local_population_index = []
uniq_k = []
for x in range(int(local_population_size/100 * biotic_components_K)):
    one = random.randint(0,biotic_components_K-1)
    while one in uniq_k:
        one = random.randint(0,biotic_components_K-1)
    uniq_k.append(one)
    local_population_index.append(one)

local_population_index.sort()

system_state = np.zeros(biotic_components_K+environment_components_N)

print(affects_w)
print(optimum_condition_u)


plt.figure(figsize=(8,8), dpi=200)
plt.title('Gaussians')
plt.xlabel('Temperature')
plt.ylabel('Biotic Effect')

for E_idx in range(environment_components_N):
    for _ in range(biotic_components_K):
        tem = []
        abundance = []
        affects = []
        abundance_sum = 0

        for temp in np.arange(0,essential_range_R, 0.001):
            tem.append(temp)
            abundance.append(((math.e) ** ((-1) * (((abs((temp)-(optimum_condition_u[E_idx][_]))) ** 2) / (2*(niche_width**2))))))
            affects.append(abundance[-1] * affects_w[0][_])
            abundance_sum += abundance[-1] * affects_w[0][_]
        plt.plot(tem,affects)

plt.show()

print(optimum_condition_u)
################## INITIAL STATE
# Abundance Init
for E_index in range(environment_components_N):
    for _ in range(biotic_components_K):
        system_state[_] = ((math.e) ** ((-1) * (((abs((En[E_index])-(optimum_condition_u[E_index][_]))) ** 2) / (2*(niche_width**2)))))
# Environment Init
for _ in range(environment_components_N):
    system_state[biotic_components_K+_] = En[_]
################## INITIAL STATE

print("System State : ", system_state)

def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change

    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    new_system_state = system_state.copy()

    for _ in range(biotic_components_K):
        if _ in local_population_index:                     # Hard coded 0 and 1 -> see notes [E1], [E1, E2], [E3][E4][E5] - scale with Es like w, u there will be another
            new_system_state[_] =  ((math.e) ** ((-1) * (((abs((system_state[biotic_components_K+0])-(optimum_condition_u[0][_]))) ** 2) / (2*(niche_width**2))))) \
                                   *\
                                   ((math.e) ** ((-1) * (((abs((system_state[biotic_components_K+1])-(optimum_condition_u[1][_]))) ** 2) / (2*(niche_width**2)))))\
                                   - system_state[_]
        else :
            new_system_state[_] =  ((math.e) ** ((-1) * (((abs((system_state[biotic_components_K+0])-(optimum_condition_u[0][_]))) ** 2) / (2*(niche_width**2))))) - system_state[_]

        #da/dt = a* - a

    for E_index in range(environment_components_N):
        biotic_force_F = 0
        for _ in range(biotic_components_K):
            biotic_force_F += (system_state[_] * affects_w[E_index][_]) # >>> contains current rate of change for alpha

        new_system_state[biotic_components_K+E_index] = (biotic_force_F)
            #dE/dt = E* + F

    return(new_system_state)

results = [[] for _ in range(biotic_components_K+environment_components_N)]

d = 0.1 # small step

x=[]
start = 0
end = 300

for step in np.arange(start, end, d):

    x.append(step)

    for _ in range(biotic_components_K+environment_components_N):
        results[_].append(system_state[_])

    k1 = d * (rates_of_change_system_state(system_state))
    k2 = d * rates_of_change_system_state(system_state + k1 * 0.5)
    k3 = d * rates_of_change_system_state(system_state + k2 * 0.5)
    k4 = d * rates_of_change_system_state(system_state + k3)

    system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)


plt.figure(figsize=(8,8), dpi=200)
plt.title('Abundance')
plt.xlabel('Time Steps')
plt.ylabel('abundance rates of change')

for _ in range(biotic_components_K):
    if _ in local_population_index:
        plt.plot(x, results[_], 'r-')
    else:
        plt.plot(x, results[_], 'b-')

plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.title('Environment Variables')
plt.xlabel('Time Steps')
plt.ylabel('Temp/PH')

plt.xlim([0, end])
plt.ylim([0, essential_range_R])

plt.plot(x, results[-2], 'k-', label = 'global')
plt.plot(x, results[-1], 'b-', label = 'local')
plt.legend()
plt.show()

