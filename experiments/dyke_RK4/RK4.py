import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import time
SAMPLE_SIZE = int(1)
SAMPLE_STEP = int(1)
RUN_ID = int(time.time())

biotic_components_K = int(100)
essential_range_R = int(100)
external_perturbation_rate_P = int(0)
time_start = int(0)
time_end = int(200)
time_step = float(1)
environment_components_N = int(1)
niche_width = int(5)
local_population_size = int(10)


################### Lotka-Volterra

def state_change_lv(state):

    new_state = state.copy() #### ELSE PYTHON PASSES IT BY REFERENCE and the loop driving this function mods its local values

    alpha = 1
    beta = 0.5
    delta = 2
    gamma = 0.5

    x = state[0]
    y = state[1]

    new_state[0] = x * (alpha - (beta * y))
    new_state[1] = -y * (gamma - (delta * x))

    return (new_state)

d = 0.001
results = [[],[]]
time = []

lv_state = np.ones(2)
lv_state[0] = 5
lv_state[1] = 5

for time_step in np.arange(0, 30, d):
    #print("start ---")
    #print("LV:",lv_state)
    lv_new = lv_state.copy()

    time.append(time_step)
    results[0].append(lv_state[0])
    results[1].append(lv_state[1])

    k1 = d * state_change_lv(lv_state)
    k2 = d * state_change_lv((lv_state + k1 * 0.5))
    k3 = d * state_change_lv((lv_state + k2 * 0.5))
    k4 = d * state_change_lv((lv_state + k3))
    lv_state += ((k1 + (2*k2) + (2*k3) + k4)/6)
    #print(lv_state)
    #print("end match start ---")


plt.figure(figsize=(8,8), dpi=200)
plt.title('LV')
plt.xlabel('time steps')
plt.ylabel('Pred + Prey')
#print("time",time)
plt.plot(time, results[0],'-')
plt.plot(time, results[1],'-')

#plt.show()



def rk4(r, h):                    #edited; no need for input f
    """ Runge-Kutta 4 method """
    k1 = h*f(r)
    k2 = h*f(r+0.5*k1)
    k3 = h*f(r+0.5*k2)
    k4 = h*f(r+k3)
    return (k1 + 2*k2 + 2*k3 + k4)/6

def f(r):
    alpha = 1.0
    beta = 0.5
    gamma = 0.5
    sigma = 2.0
    x, y = r[0], r[1]
    fxd = x*(alpha - beta*y)
    fyd = -y*(gamma - sigma*x)
    return np.array([fxd, fyd], float)

h=0.001                               #edited
tpoints = np.arange(0, 30, h)         #edited
xpoints, ypoints  = [], []
r = np.array([2, 2], float)
for t in tpoints:
    xpoints.append(r[0])          #edited
    ypoints.append(r[1])          #edited
    r += rk4(r, h)             #edited; no need for input f

plt.plot(tpoints, xpoints)
plt.plot(tpoints, ypoints)
plt.xlabel("Time")
plt.ylabel("Population")
plt.title("Lotka-Volterra Model - no T")
plt.savefig("Lotka_Volterra.png")
#plt.show()

################### Lotka-Volterra


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

system_state = np.zeros(biotic_components_K+1)

print(affects_w)
print(optimum_condition_u)

plt.figure(figsize=(8,8), dpi=200)
plt.title('Gaussians')
plt.xlabel('x')
plt.ylabel('y')


for _ in range(biotic_components_K):
    tem = []
    abundance = []
    affects = []
    abundance_sum = 0

    for temp in np.arange(0,essential_range_R, 0.001):
        tem.append(temp)
        abundance.append(((math.e) ** ((-1) * (((abs((temp)-(optimum_condition_u[0][_]))) ** 2) / (2*(5**2))))))
        affects.append(abundance[-1] * affects_w[0][_])
        abundance_sum += abundance[-1] * affects_w[0][_]
    plt.plot(tem,affects)

plt.show()

################## INITIAL STATE
E = 50
system_state[-1] = E
for _ in range(biotic_components_K):
    system_state[_] = ((math.e) ** ((-1) * (((abs((E)-(optimum_condition_u[0][_]))) ** 2) / (2*(5**2)))))
################## INITIAL STATE

print("System State : ", system_state)

def rates_of_change_system_state(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change

    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    new_system_state = system_state.copy()

    for _ in range(biotic_components_K):
        new_system_state[_] =  ((math.e) ** ((-1) * (((abs((system_state[-1])-(optimum_condition_u[0][_]))) ** 2) / (2*(5**2))))) - system_state[_]
        #da/dt = a* - a

    biotic_force_F = 0
    for _ in range(biotic_components_K):
        biotic_force_F += (system_state[_] * affects_w[0][_]) # >>> contains current rate of change for alpha

    new_system_state[-1] = (biotic_force_F)
        #dE/dt = E* + F

    return(new_system_state)

results = [[] for _ in range(biotic_components_K+1)]

d = 0.1 # small step

x=[]
start = 0
end = 400

biotic = []

for step in np.arange(start, end, d):

    x.append(step)

    for _ in range(biotic_components_K+1):
        results[_].append(system_state[_])

    force = 0
    for _ in range(biotic_components_K):
        force += (system_state[_] * affects_w[0][_])
    biotic.append(force)

    k1 = d * (rates_of_change_system_state(system_state))
    k2 = d * rates_of_change_system_state(system_state + k1 * 0.5)
    k3 = d * rates_of_change_system_state(system_state + k2 * 0.5)
    k4 = d * rates_of_change_system_state(system_state + k3)

    system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)


plt.figure(figsize=(8,8), dpi=200)
plt.title('Alpha Rates')
plt.xlabel('x')
plt.ylabel('alpha rates')

for _ in range(biotic_components_K):
    plt.plot(x, results[_], '-')

plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.title('Environment Variable')
plt.xlabel('x')
plt.ylabel('E')

#plt.xlim([0, 100])
#plt.ylim([-100, 100])

for _ in range(len(x)):
    plt.plot(x, results[-1], '-')
plt.show()

plt.figure(figsize=(8,8), dpi=200)
plt.title('Biotic Force')
plt.xlabel('x')
plt.ylabel('Biotic Force')
plt.plot(x, biotic, '-')

plt.show()


