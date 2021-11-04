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
    SAMPLE_SIZE = s['SAMPLE_SIZE']
    SAMPLE_STEP = s['SAMPLE_STEP']
    RUN_ID = s['RUN_ID']

    biotic_components_K = s['biotic_components_K']
    essential_range_R = s['essential_range_R']
    external_perturbation_rate_P = s['external_perturbation_rate_P']
    time_start = s['time_start']
    time_end = s['time_end']
    time_step = s['time_step']
    environment_components_N = s['environment_components_N']
    truncated_gaussian_ROUND = s['truncated_gaussian_ROUND']
    niche_width = s['niche_width']
    local_population_size = s['local_population_size']
    affects_w = s['affects_w']
    optimum_condition_u = s['optimum_condition_u']
    biotic_force_F = s['biotic_force_F']
    local_population_index = s['local_population_index']

    global_start_temp = s['global_start_temp']
    local_start_temp = s['local_start_temp']

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']

finally:
    s.close()

S_STEP = SAMPLE_STEP
K = biotic_components_K
R = essential_range_R
P = external_perturbation_rate_P
start = time_start
end = time_end
step = time_step
N = environment_components_N
E = [global_start_temp, local_start_temp]
F = biotic_force_F
print(E)

ROUND = truncated_gaussian_ROUND

OEn = niche_width
OE = [OEn for _ in range(K)]

w = affects_w
u = optimum_condition_u

#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]

rNumberAlive = [[] for x in range(K)]

# After three steps - the equilibrium abundance stagnates - hence loop below till 5
rEquilibrium_Abundance = [[] for x in range(K)]

rStart_Alives = [[] for x in range(K)]

#for _ in range(K):
#    rNumberAlive[_].append(0)

alpha = [[] for _ in range(K)]

# Fix applied - correctly generating alphas for global and local
for _ in range(K):
    al = []
    for ai in range(N):
        al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

        new_alpha = 0
        if _ in local_population_index:
            new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
        else:
            new_alpha = al[0]       # Else take the first one as Eg

        alpha[_].append(new_alpha)

        rAx[_].append(new_alpha)
        rAxR[_].append(new_alpha * R)
        rEquilibrium_Abundance[_].append(new_alpha)

        if(new_alpha > 0):
            rNumberAlive[_].append(1)
            rStart_Alives[_].append(1)
        else:
            rNumberAlive[_].append(0)
            rStart_Alives[_].append(0)

rF = [[] for _ in range(N)]         #Biotic Force Values
rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

initfSUM = [0 for _ in range(N)]
for _ in range(K):
    initfSUM[0] = initfSUM[0] + (alpha[_][-1] * w[0][_])

rF[0].append(initfSUM[0])

for _ in range(K):
    if( _ in local_population_index):
        initfSUM[1] = initfSUM[1] + (alpha[_][-1] * w[1][_])

rF[1].append(initfSUM[1])


for _ in range(N):
    rE[_].append(E[_])                 #Input the Start Temperatures

#Tracks time (steps accumulation)
time = []

biotic_force = [[] for _ in range(K)]
temperatures = []


def update(step):
    global F, P, E, Et, rF, rP, rE, rEt, u, w, N

    fSUM = [0 for _ in range(N)]
    alpha_time_scale = 0.7
    temperature_time_scale = 0.2

    for _ in range(K):
        al = []
        for ei in range(N):
            al.append(round((math.e) ** ((-1) * (((abs((E[ei])-u[ei][_])) ** 2) / (2*(OE[_]**2)))), ROUND))

        new_alpha = 0
        if( _ in local_population_index):
            new_alpha = np.prod(al)
        else:
            new_alpha = al[0]

        k1 = new_alpha
        k2 = k1 + (k1 * step/2)
        k3 = k1 + (k2 * step/2)
        k4 = k1 + (k3 * step)
        yt = alpha[_][-1] + (((k1 + (2*k2) + (2*k3) + k4)/6) - alpha[_][-1]) * step

        alpha[_].append(yt)

        rAx[_].append(yt)
        rAxR[_].append(yt * R)

        if(yt > 0):
            rNumberAlive[_].append(1)
        else:
            rNumberAlive[_].append(0)

    for _ in range(K):
        fSUM[0] = fSUM[0] + (alpha[_][-1] * w[0][_])

    F[0] = fSUM[0]
    newE = E[0] + (((0 + F[0]) * 1))

    # LOCKING E No Change to test Equilibrium
    E[0] = E[0] #+ ((newE-E[0]) * temperature_time_scale)
    rF[0].append(F[0])
    rE[0].append(E[0])

    # ============ END EG ==============

    for _ in range(K):
        if( _ in local_population_index):
            fSUM[1] = fSUM[1] + (alpha[_][-1] * w[1][_])

    F[1] = fSUM[1]
    newE = E[1] + (((0 + F[1]) * 1))

    # LOCKING E No Change to test Equilibrium
    E[1] = E[1] #+ ((newE-E[1]) * temperature_time_scale)
    rF[1].append(F[1])
    rE[1].append(E[1])

    # ============ END EL ==============

def results_shelve():

    s = shelve.open(shelve_file)

    try:
        s['rEquilibrium_Abundance'] = rEquilibrium_Abundance
        s['rStart_Alives'] = rStart_Alives

    finally:
        s.close()

if __name__ == '__main__':

    #time.append(0)
    post_init_start = start + step

    for xtime in np.arange (post_init_start, post_init_start+5, step):
        update(step)
        time.append(xtime)

    results_shelve()












