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
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
plt.rcParams["font.family"] = "Times New Roman"

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

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
ENV_START=[5]
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
    plt.title('Biotic Force of 100 species with alive threshold : ' + str(ALIVE_THRESHOLD), fontsize=30)
    plt.xlabel('Temperature', fontsize=20)
    plt.ylabel('Biotic Force', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(SPECIES_K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
    plt.legend(loc='upper right')
    plt.show()

def fYa(Xe, Ni, u):
    return (((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2))))))

def fXe(Ya, Ni, u):
    return (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))) + u

def plot_inverse():

    Xe = 60 # spot on x axis - temperature
    Ya = 0  # abundance
    Ni = 5 # niche
    u = 50 # ideal growing temperature

    #Ya = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))
    #Xe = (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))) + u
    #Xe = np.linspace(-50,150,1000)


    f2 = np.vectorize(fYa)
    x = np.arange(-5, 105, 0.001)
    plt.plot(x, f2(x, Ni, u))
    plt.show()

    f2 = np.vectorize(fXe)
    x = np.arange(0.00001, 1, 0.0001)
    plt.plot(f2(x, Ni, u), x)
    plt.show()


    print(fYa(58.97061288997051, Ni, u))
    print(fXe(1, Ni, u))
    print(fXe(0.2, Ni, u))

#plot_inverse()


def fYaI(Xe, Ni, u, T):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if(abundance <= T):
        abundance = 0

    return(abundance)

def fXe(Ya, Ni, u):
    return (u + (math.sqrt(((-1 * math.log(Ya,math.e)) * (2*(Ni**2))))))

def fXe_negative(Ya, Ni, u):
    return (u - (math.sqrt(((-1 * math.log(Ya,math.e)) * (2*(Ni**2))))))

def fYaIx(Xe, Ni, u, NRange):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if((Xe >= u + NRange) or (Xe <= u - NRange)):
        abundance = 0

    return(abundance)

def plot_inverse_case():

    Ni = 5 # niche
    u = 50 # ideal growing temperature
    ST = 0.2

    #print(fYaI(45, Ni, u, ST))

    #f2 = np.vectorize(fYaI)
    #x = np.arange(-5, 105, 0.001)
    #plt.plot(x, f2(x, Ni, u, ST))
    #plt.show()

    alpha = []
    temps = []

    for x in np.arange(-5, 105, 0.001):
        temps.append(x)
        alpha.append(fYaI(x, Ni, u, ST))

    fig, ax = plt.subplots(figsize=(20,20), dpi=300)
    ax.set_xlabel('Temperature', fontsize=30)
    ax.set_ylabel('Abundance', fontsize=30)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1, labelsize=20)
    ax.tick_params(which='major', length=8, labelsize=20)
    ax.tick_params(which='minor', length=6, labelsize=20)
    ax.axhline(y=0.2, color='r', linestyle='-', label='survival threshold')
    ax.plot(temps,alpha)
    ax.legend(prop={'size': 30})
    fig.tight_layout()
    fig.savefig("ST Model")
    fig.show()





    Ni = 20

    alpha = []
    temps = []

    for x in np.arange(-5, 105, 0.001):
        temps.append(x)
        alpha.append(fYaI(x, Ni, u, ST))



    fig, ax = plt.subplots(figsize=(20,20), dpi=300)
    ax.set_xlabel('Temperature', fontsize=30)
    ax.set_ylabel('Abundance', fontsize=30)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1, labelsize=20)
    ax.tick_params(which='major', length=8, labelsize=20)
    ax.tick_params(which='minor', length=6, labelsize=20)

    ax.axhline(y=0.2, color='r', linestyle='--', label='survival threshold 0.2')
    ax.plot(temps,alpha,'--', color="#FF7F0E", label='IGAF Model 0.2')

    #plt.tight_layout()
    #plt.show()

    #print(fXe(0.2, 5, u))

    NRange = (fXe(0.2, 5, u) - u)

    #print("NRANGE")
    #print(NRange)

    alpha = []
    temps = []


    for x in np.arange(-5, 105, 0.001):
        temps.append(x)
        alpha.append(fYaIx(x, Ni, u, NRange))

    ax.axhline(y=fYaI(fXe(0.2, 5, u), 20, u, 0), color='r', linestyle='-', label='survival threshold 0.9')
    #plt.annotate('truncation : '+str(fYaI(fXe(0.2, 5, u), 20, u, 0)), xy =(0.5, 1), xytext =(0.5, 0.5))
    print("Truncation : " + str(fYaI(fXe(0.2, 5, u), 20, u, 0)))

    ax.plot(temps,alpha, label='IGAF Model 0.9')

    handles, labels = plt.gca().get_legend_handles_labels()

    #specify order of items in legend
    order = [1,0]

    #add legend to plot
    ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
    ax.legend(prop={'size': 30}, loc="center right", bbox_to_anchor=(1, 0.7))
    fig.tight_layout()
    fig.savefig("ST Model")
    fig.show()



    #print("Truncation : " , str(fYaI(fXe(0.2, 5, u), 20, u, 0)))

    #print(fXe(0.5, 20, u))
    #print(fXe_negative(0.5, 20, u))

    #print("Verification : ")

    #print(fYaI(fXe(0.5, 20, u), 20, u, 0))
    #print(fYaI(fXe_negative(0.5, 20, u), 20, u, 0))



truncation = 0.2


def plot_temps():

    temperatures = []
    step = 0.01

    for x in np.arange (-25, RANGE_R+25, step):
        temperatures.append(x)

    abundance_for_st = [[] for _ in range(SPECIES_K)]
    step = 0.01

    for y in range(SPECIES_K):
        for x in np.arange (-25, RANGE_R+25, step):
            aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
            if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                abundance_for_st[y].append(0)
            else:
                abundance_for_st[y].append(aliveness)

    abundance_for_st_x = [[] for _ in range(SPECIES_K)]

    for y in range(SPECIES_K):
        for x in np.arange (-25, RANGE_R+25, step):
            NRange = (fXe(0.2, 5, mu[0][y]) - mu[0][y])
            aliveness = fYaIx(x, 10, mu[0][y], NRange)
            abundance_for_st_x[y].append(aliveness)


    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))

    #fig.set_size_inches(3, 1.5)
    for _ in range(SPECIES_K-80):
        ax1.plot(temperatures,abundance_for_st[_])
    #ax1.set_title('ST Model', fontsize=20)

    ax1.set_xlabel('Temperature', fontsize=30)
    ax1.set_ylabel('Abundance', fontsize=30)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(which='both', width=1)
    ax1.tick_params(which='major', length=8)
    ax1.tick_params(which='minor', length=6)

    ax1.set_xlabel('Temperature', fontsize=20)
    ax1.set_ylabel('Abundance', fontsize=20)
    for _ in range(SPECIES_K-80):
        ax2.plot(temperatures,abundance_for_st_x[_])
    #ax2.set_title('NW Model', fontsize=20)
    ax2.set_xlabel('Temperature', fontsize=20)
    ax2.set_ylabel('Abundance', fontsize=20)

    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(which='both', width=1)
    ax2.tick_params(which='major', length=8)

    fig.tight_layout()
    fig.show()


def plot_alphas_truncated():


    temperatures = []
    biotic_force = [[] for _ in range(SPECIES_K)]
    step = 0.01

    for x in np.arange (-25, RANGE_R+25, step):
        temperatures.append(x)

    for y in range(SPECIES_K):
        for x in np.arange (-25, RANGE_R+25, step):
            aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
            if(aliveness <= truncation and aliveness >= (-1 * truncation)):
                biotic_force[y].append(0)
            else:
                biotic_force[y].append(aliveness * omega[0][y])



    fig, ax = plt.subplots(figsize=(20,20), dpi=300)
    ax.set_xlabel('Temperature', fontsize=30)
    ax.set_ylabel('Biotic Force', fontsize=30)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1, labelsize=20)
    ax.tick_params(which='major', length=8, labelsize=20)
    ax.tick_params(which='minor', length=6, labelsize=20)
    for _ in range(SPECIES_K):
        ax.plot(temperatures,biotic_force[_])

    ax.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4, label='Combined Biotic Force')
    ax.legend(prop={'size': 30})
    fig.tight_layout()
    fig.savefig("ST Model")
    fig.show()


    temperatures = []
    step = 0.01

    for x in np.arange (-25, RANGE_R+25, step):
        temperatures.append(x)


    abundance_for_st_x = [[] for _ in range(SPECIES_K)]

    for y in range(SPECIES_K):
        for x in np.arange (-25, RANGE_R+25, step):
            NRange = (fXe(0.2, 5, mu[0][y]) - mu[0][y])
            aliveness = fYaIx(x, 10, mu[0][y], NRange)
            abundance_for_st_x[y].append(aliveness * omega[0][y])


    fig, ax = plt.subplots(figsize=(20,20), dpi=300)
    ax.set_xlabel('Temperature', fontsize=30)
    ax.set_ylabel('Biotic Force', fontsize=30)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1, labelsize=20)
    ax.tick_params(which='major', length=8, labelsize=20)
    ax.tick_params(which='minor', length=6, labelsize=20)
    for _ in range(SPECIES_K):
        ax.plot(temperatures,abundance_for_st_x[_])
    ax.plot(temperatures,np.sum((np.array(abundance_for_st_x, dtype=float)), axis=0), lw=4, label='Combined Biotic Force')
    ax.legend(prop={'size': 30})
    fig.tight_layout()
    fig.savefig("ST Model")
    fig.show()



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
    alive_thresh = 0.2
    for each_temp in np.arange(0,100,0.01):
        temp_t.append(each_temp)
        result = (math.e) ** ((-1) * (((abs(each_temp-ideal_temp)) ** 2) / (2*(NICHE**2))))
        if (result > alive_thresh):
            gaus_t.append(result)
        else:
            gaus_t.append(0)

    #plt.axhline(y=ALIVE_THRESHOLD, color='g', linestyle='--')
    #plt.plot(temp_t,gaus,_t 'b-',label = 'The gaussian distribution')
    #plt.show()


    fig, (ax1, ax2) = plt.subplots(1, 2, dpi=300, figsize=(30,10))
    fig.suptitle('The survival threshold',fontsize=30)
    #fig.set_size_inches(3, 1.5)
    ax1.plot(temp, gaus)
    ax1.set_title('The Original Model', fontsize=20)
    ax1.set_xlabel('Temperature', fontsize=20)
    ax1.set_ylabel('Abundance', fontsize=20)
    ax2.plot(temp_t, gaus_t)
    ax2.set_title('Survival Threshold of 0.2', fontsize=20)
    ax2.set_xlabel('Temperature', fontsize=20)
    ax2.set_ylabel('Abundance', fontsize=20)

    ax2.hlines(y=0.2, xmin=0, xmax=100, linewidth=2, color='r', label='survival threshold')
    ax2.legend(prop={'size': 20})
    fig.show()


def rates_of_change_system_state_new_niche(system_state):

    # Environment Vars Change >>> Abundance >>> Biotic Force Changes >>> Environment Vars Change\
    # Alphas_IN determine E_OUT via biotic Force
    # E_IN determine Alphas_OUT via Gaussian

    rate_of_change = system_state.copy()

    Eg = system_state[SPECIES_K+0]


    for s_i in range(SPECIES_K):
        NRange = (fXe(0.2, 5, mu[0][s_i]) - mu[0][s_i])
        a_star = fYaIx(Eg, 10, mu[0][s_i], NRange)
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

    for sim in range(0, 300):

        #print(sim)
        omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        print(sim)
        print(omega)
        print(mu)
        #omega = [[-0.7606309166232983, -0.30656459261181834, 0.05129001309671177, 0.8222774249661196, -0.542971708146222, -0.6568450067100073, -0.2496095511396157, -0.2668532256929059, 0.8881170192804011, 0.909508485820103, -0.20755931188143273, -0.49525968173430934, 0.7146094503469982, 0.858618177871677, -0.2756183587753862, 0.11248630329649356, -0.598689670721779, 0.7053711716887376, 0.28139746988367964, -0.8512075764249425, -0.7516160270573229, -0.7500594499640028, -0.9903407710753105, -0.4071967617633194, 0.9326976323171305, 0.6696594573582442, 0.5615278795832817, 0.0018042188896614242, 0.14364585275926034, -0.1475146556195548, -0.9818272182803782, -0.9606105204746622, 0.5521104337156226, -0.985524212027862, -0.4501535347912844, 0.7347080500639207, -0.48687751725120787, -0.25010603033145484, -0.7260614787652673, 0.18522993202870142, 0.2613221089767306, -0.9597428598835354, 0.49120311384107884, 0.9133916790827177, -0.5270950170393662, -0.7091124684269821, -0.06438304709440201, -0.09317320894962888, -0.9492938163564253, -0.49360237913426785, 0.698926431716375, 0.5901097772618877, -0.5179239893046841, 0.3854529991583602, 0.13240463789438417, 0.5081645866440943, 0.03433347412512022, 0.5549361572895568, -0.981220705141056, 0.7476002079344164, 0.978180008981473, -0.2037893253909986, 0.9333536806382772, 0.6266215863870039, 0.9135455527390071, -0.16301660472149537, 0.0990767640653456, -0.33110844720401644, 0.9438425879580892, 0.9274388294755305, 0.29873647062917463, 0.5817388966871686, -0.0043540073620593756, 0.49597487507407156, 0.3443365526042439, 0.1430984387907408, -0.7502893720684516, 0.36144670690924396, 0.33977413082691155, -0.5918106904314755, 0.7921032441676021, 0.04712772558509215, 0.6501512170529951, 0.3983894753380126, 0.6520404217280338, 0.018104903738153988, 0.6881594438436913, -0.05544429899984227, 0.3951868524550115, 0.1254238110056256, 0.21741232043919156, 0.894811789615579, 0.25042489666576095, -0.27455410732584595, -0.495371342773433, -0.7295573673764788, -0.3164041835355933, -0.8861918776774904, 0.5161605204331021, -0.31997066940613283]]
        #mu = [[27.47072653065391, 72.23852775955119, 16.62261751218004, 28.609485834284744, 73.01993304058304, 30.24013979458361, 57.23888786373822, 51.75784515865691, 86.23481442996228, 79.38739736743509, 26.76485401953976, 79.86784620321485, 87.19163136649745, 37.38157055724236, 60.65223516627422, 90.99013728695155, 67.94046108980127, 0.09934338840221768, 27.46472652440971, 4.442104570899874, 74.4454301381862, 43.18903024166305, 65.41425056367859, 21.772941738257035, 79.0146426820896, 73.52028427485699, 73.69158268909696, 51.65306846474306, 34.24965413161989, 44.852275128847474, 41.14957583626929, 7.773295206297581, 91.56228670579549, 62.177459233333146, 16.563610645324257, 69.36913941271399, 94.92474590221734, 10.4565993024903, 71.34167763108775, 81.18561539523029, 10.659502267218034, 59.48452232340781, 1.434446943775869, 25.935072755774545, 91.80653805141665, 77.433236458811, 16.245524578050208, 17.473934910891032, 44.96212363168001, 84.61360175038361, 45.43519383163412, 95.25995748291835, 47.592461034809794, 14.306512034645813, 68.54317884776653, 35.06874402868364, 49.9896801520091, 78.45256083336089, 80.76923451949591, 39.67799061972166, 13.037849998417927, 73.96845064664885, 35.64790567442764, 92.3768308177704, 36.548940847278125, 10.856549812549531, 92.63888936119639, 86.3373773278121, 8.091107112224648, 92.94483288818788, 62.47513902342553, 75.93483106082417, 82.00606718997655, 33.098261036676426, 22.15594427323587, 90.88174468441669, 72.72053816890872, 46.53506895978241, 74.41126995089272, 45.19184314310557, 12.482157877347388, 77.55972523656204, 91.84338411933805, 9.307824392803244, 85.37223780735305, 87.572454037445, 98.97077903829003, 15.450437287172136, 85.90015532127354, 31.68767781601015, 5.6201832599341195, 32.18538022870965, 14.539033055396033, 99.89448344574492, 84.04861338729293, 42.76573523713013, 70.01653986983355, 24.304995167303634, 96.22304536772086, 90.72238883117949]]

        omega = [[-0.05098490038776626, 0.17455678984045542, 0.26756998137923493, 0.13502937585932706, -0.974690530107956, -0.19642938829181156, -0.9709124945974248, 0.6821713333948187, 0.7804520470172993, 0.5046923274805883, 0.3426430338956836, -0.2593501991035976, -0.3447428278296396, -0.37651138829194997, 0.6802538792552002, -0.18511885004273299, 0.8535909151584309, 0.6838330334797145, -0.19557245241836818, -0.6322141128030832, -0.4097378512002394, -0.717058203715834, -0.36148725471643495, -0.07063197041983749, 0.8425553546433817, -0.4715024535051524, 0.83869264996091, 0.9258465879612079, 0.008790445100292477, -0.15642600873813772, 0.16127030834314104, 0.9066052467583974, 0.5342354826558324, 0.6521759006378407, -0.01466324405255115, -0.8038781010641236, -0.9697471909770907, 0.5547883439812757, -0.0848743423716618, -0.7796248296964636, -0.7920376588136944, 0.1627834086172537, -0.002600307858052764, -0.7349295989768028, 0.2105139868276802, 0.14291113890004703, 0.3878480976033305, -0.92367878949198, 0.6828237004076034, -0.3424080628942132, 0.44393466393905756, -0.3196425123555642, 0.605554475428342, 0.6881694602967869, 0.1464848525921052, 0.2713213263421108, 0.15366524241870527, 0.0037674957663733633, -0.9836574803467097, -0.42930529771499604, -0.07912277397166267, -0.7849124139504611, -0.3017105728768146, -0.16973003839555334, 0.8685003713879662, -0.9902985976631122, 0.6227883571339965, -0.8593059180057268, 0.14123251697249373, 0.6739848383876017, 0.18985660959767836, 0.09734930302979272, -0.14517832364294692, 0.6628792172342415, -0.08047355752221952, 0.6164234151789485, 0.38911747267993535, -0.42595520632983064, 0.63141871079892, -0.9366311127428604, -0.6904017980332149, 0.16877390370683454, 0.23658075007687018, 0.00904017946067981, -0.17365702264196958, 0.6877721002667452, 0.42097814810810186, -0.32770664077129674, -0.9232255253715389, -0.5590063526681543, 0.8606653417089609, -0.26265996326004437, 0.22611150366084698, 0.23732649361435354, 0.956834070966317, 0.8437773250552523, -0.13076144352119168, -0.7973834641142894, 0.9404315112684978, 0.29853011648480243]]
        mu = [[82.24089937860853, 54.17316120627653, 77.56536741045093, 84.03752589868849, 82.92046645302341, 17.0738947577205, 26.71564056208009, 14.735337283852001, 81.58139855104717, 76.11620527210489, 58.704930913180895, 50.42506840059533, 24.95007136098343, 26.693011789547995, 50.68267162856876, 79.93846610081862, 6.874905812788512, 67.06934373821062, 53.31220316850078, 18.198458427491037, 18.119553178164548, 78.04238990575159, 20.72357884480137, 25.76454666552509, 16.722233500414585, 45.93151429091162, 70.72676803957418, 52.006021020376245, 95.24806799376067, 92.46831068854327, 84.41157337941382, 53.96643491074653, 87.45211052792946, 25.527559251881993, 93.15289961815463, 75.51351585955287, 38.848090589588836, 54.98207551343427, 15.0554359935822, 71.26805172347794, 17.429465903823004, 31.795963719656072, 70.5105418355675, 68.47557996333803, 53.03458833553273, 43.74044464075689, 1.6812522312451894, 34.471577809802525, 34.542255418696186, 95.03267098786657, 3.337357483175385, 22.680020626028707, 17.327624404216515, 45.22396313895687, 89.33969098233804, 79.51396640786224, 66.24690658491626, 15.371410425560384, 69.47813611883967, 1.5984414704928174, 64.73139798907268, 35.052389676064834, 39.140226107076884, 79.14548487675731, 92.34300547767702, 88.05314557465054, 88.84999283654497, 89.06458188064067, 91.6209260057628, 73.33000191344694, 74.51763603497695, 63.72004151791714, 10.240002620952271, 4.993476058516045, 84.38848271297779, 55.88896053324516, 49.79736985609149, 77.32066757789754, 99.29231342511862, 37.62061967370821, 4.460698337070945, 63.89066777413285, 62.238558706921154, 75.16223865133166, 51.50727476915894, 82.47584554107272, 57.98590168812551, 12.335985714962927, 9.725501648974522, 20.735157859452624, 31.841015794832895, 55.79903726187043, 65.41899056394642, 33.93865828850213, 48.11869704913487, 93.86064794946483, 81.24043115242094, 55.58538285232855, 50.810533456768155, 33.04038229916346]]


        # For Write Up
        # NW being better at capturing edge
        # # 143
        #omega = [[-0.49903771532678776, 0.40913226401524194, 0.10387566381122637, 0.03960774949455503, 0.926379746461423, -0.4399989329366092, -0.9965802079712764, -0.577142656525907, 0.9446806538291668, 0.300318792432857, 0.9873286288874199, 0.12307041127091889, 0.7998548140198871, -0.556711082702406, -0.5942002992449338, -0.056234075638191605, 0.8337771359475838, -0.0907109801903121, 0.47049127967275406, -0.046227242243688904, -0.8605054847221216, -0.3842889218994212, -0.3683921061809805, -0.9729985920882773, -0.9438039325735446, 0.00727188484057284, 0.8069639201530387, 0.4218790948715676, 0.6443782099291055, 0.2883768999641094, -0.46473137314355784, -0.010740539664183535, -0.7704380148141969, -0.6018760114183308, -0.32168749462177026, 0.5399650444353303, -0.7434431547828924, -0.5483432120086935, 0.44231304443527053, 0.2346416431997982, 0.6445584960022823, -0.38850026746471955, -0.21520672109792205, -0.000159761403309977, -0.7717844655957851, -0.5432465491906262, -0.8216563204218548, -0.16521328372181765, -0.1924323235975547, 0.14816184650979958, -0.049414783398535, -0.6969891571340863, -0.07894640422914367, -0.6927921917657891, -0.7644838072444173, -0.2962900801615409, -0.9156125987227017, -0.47668179135512956, -0.05241808415345095, 0.6987600874085629, -0.8396597479139323, 0.3929609652516166, -0.41998294959820237, -0.07367050490862392, -0.27252434224908795, -0.6326276754574751, 0.8552465371059552, 0.518184467575753, -0.5018766937728643, 0.5166510779700213, -0.8839856153268912, -0.41471162810165474, -0.43298182468426116, -0.4115765361554149, 0.6342381333326106, -0.3399305411060054, 0.9797117417031356, 0.30910059557010006, 0.8707829951149579, -0.48855507359857997, -0.4676215272346944, 0.8990930775633987, -0.06551544648249874, 0.40342130250347275, -0.3818114415053526, -0.8269758135322605, -0.4511904620486491, -0.0816814732036486, -0.35332893451440706, -0.6004113686412205, -0.7490005469131635, 0.5474328441528709, -0.21445216595468453, -0.8590098267005353, -0.530290095729983, 0.940681723065034, 0.08991020650707071, -0.37477113592727296, 0.3804119044854015, 0.782847654180816]]
        #mu = [[29.594855868804782, 87.45486635152619, 32.1588020735042, 38.058134180395456, 13.176026301495403, 61.42070551419646, 64.20900796697863, 84.65959410638598, 75.10122106190782, 14.314327791284231, 81.16639885571439, 37.67362979309158, 56.99093635436483, 60.265094921114304, 49.89864307565098, 98.63918627919712, 20.054350442705427, 23.49214301845346, 11.632683955722879, 27.831315846645076, 38.38154693675517, 50.17795033821499, 94.66199385277238, 61.5106547547274, 95.44853579180685, 6.133363050705332, 27.287722565461458, 22.934913927049827, 5.553891961338897, 93.22402599266583, 76.07968848598048, 73.58520623938657, 66.72479722400485, 5.826723695883718, 79.78116237319784, 26.22123680817947, 87.06966146244784, 58.405161394901896, 71.66762494264931, 49.80642558541454, 23.901096602512318, 48.703114631131385, 89.76437979024507, 68.02235036135264, 63.15784845164877, 55.8905242198511, 34.54634924117771, 94.51951780085717, 93.54099718271584, 75.58073357118604, 38.44003969987724, 40.65908770427917, 26.062687422600717, 1.794391319548394, 85.61805870366722, 91.04171607181144, 60.53513043817564, 22.7227508760978, 71.14371951796576, 41.8698480656971, 25.64608854900088, 97.52802698022565, 71.86158147821041, 33.89182765439538, 96.37028759209882, 34.830247115922774, 5.408312764490653, 56.88324488766654, 9.889238160950576, 27.32151696750532, 97.55860396857233, 43.155670200701394, 5.085733098820144, 27.62205300337679, 53.49969270942497, 28.41784711264245, 98.71323433364896, 60.6610286787993, 67.5886345684845, 71.20717222579513, 58.24215386992652, 85.92517484011614, 48.09859193273523, 18.761945393763234, 43.567086240332095, 36.9604419893391, 16.122288327894186, 71.81253549059876, 91.30855956365946, 55.4365257139535, 28.899033132169873, 64.35534926932618, 76.1728579362988, 97.77438285155698, 88.32092829879569, 31.647453593576714, 23.37186519599873, 7.872480422149419, 69.19095153335229, 85.73153107855258]]
        # #
        # # 161
        # omega = [[0.7118125706605243, -0.5404464058892815, -0.012156390161768948, 0.4926345587000993, -0.07447508455125695, -0.3399152240993284, 0.6413989459367888, -0.010173698601906978, -0.8160155411697856, 0.9112397896263214, 0.8215206925161356, -0.1451703224845733, 0.32378108752644486, -0.9010282324251384, 0.8646403455451124, 0.8686730202415258, -0.2447201526655356, 0.9968744345941314, 0.5883437218820455, 0.3810518045454159, 0.09403042188893651, -0.7598987608166436, -0.48607532764343264, -0.6874871448784976, 0.11406269045650874, -0.511891568251257, -0.21930206826984033, -0.7125430273574471, -0.4547403343769465, -0.9865884316365827, -0.18241947379959011, -0.16873260148652713, 0.847105149968697, -0.9366868336531089, 0.6613547770663613, 0.5908754418685866, 0.6105995120743841, 0.2656674725637582, 0.09380050832002063, 0.27217851696383377, -0.7759611392644972, 0.9144537525091709, -0.4113271723117067, -0.13153042382811586, -0.511480504531304, -0.284911841687719, 0.14253663931228133, 0.8653051986337121, -0.7004348625905896, 0.3382305228154563, -0.8908469769309888, 0.8110257176965168, -0.35372380832701666, -0.6406150892185913, -0.1869035244535291, -0.17555321538504742, -0.058824667186782253, 0.38931868569293826, -0.8158461528367986, 0.9983037263422843, -0.2569070167797016, -0.5506400917043468, 0.4822722821021592, -0.3034860678521134, 0.2223923804717991, 0.11153747976996686, 0.23884558636959152, -0.6787067730965926, 0.9515701482687837, 0.5615611750003178, -0.3627966227467212, -0.15617299400098394, 0.9570533677148261, 0.04749584558872222, -0.9064972152821233, -0.3258760195754111, -0.1978340265258023, -0.427519680601939, -0.8204753205387354, 0.22476708355331043, 0.9710036620873475, -0.6262202912627153, 0.7623206197981396, -0.8326789276529518, -0.11390161912581598, 0.018845097675205258, 0.8540972027022644, -0.3457779484126775, -0.3252266326435236, -0.8179527684118062, 0.5413507169150795, 0.4469286689365648, 0.4569970970822661, 0.4974796340388319, -0.48932184130465495, 0.7387411564661688, -0.6048295327286959, -0.5364659580377062, -0.9912754859969519, -0.21432479230862933]]
        # mu = [[77.27997357460593, 89.25908700482773, 92.02542812283625, 5.497061870977616, 97.63867487342691, 46.33259371773435, 59.88633919068431, 45.37147186044059, 72.34022527075146, 54.22872278671612, 97.83958904884679, 28.22274574056428, 95.69345221099933, 71.31154974899705, 18.631444982380007, 89.74969550166134, 48.12913722942183, 9.078384149681495, 20.389455416882818, 98.30371342535503, 19.12793397277707, 70.62469765102463, 53.49862161458191, 53.827750903838776, 44.27345519315667, 45.60133105773802, 40.59289176856801, 29.03394019778616, 69.54454418973873, 41.184404710627746, 12.676316958823563, 11.02729966434609, 14.188602409223428, 31.42007603951784, 80.07507523172315, 9.275196148519472, 89.44082752066421, 0.7177986349034771, 19.46563266675787, 27.466388306387614, 92.41208853716158, 99.31475624851404, 5.241795339494848, 52.20491054573102, 73.30576391525112, 27.28341856476739, 63.08271257999858, 23.857187598021678, 3.8084653359235054, 57.442891700515034, 93.06572112208814, 63.030282145916026, 28.55349691712784, 29.66896870345861, 37.290997342963635, 71.39842258310102, 2.145856418121439, 5.893402441139117, 61.44788387927035, 56.10152596190549, 51.880862112790346, 11.569522536930366, 91.84062216632667, 41.34896555555086, 20.151447238125186, 86.36825470815705, 61.63274014886837, 92.50556725192213, 28.116816077646757, 76.89605071559603, 53.13976026800338, 13.841178770746343, 46.311078378134575, 85.58150499878336, 53.28874993554172, 90.73440454109434, 72.38741824023309, 73.65837209991014, 61.116897759432156, 77.55293376398916, 46.318739532572565, 16.81476707944678, 58.88560348543032, 54.487213797933485, 31.15873645732069, 70.13773680796017, 93.56815539320851, 0.5464990027113092, 77.27767921305623, 5.780980272898939, 61.11828591508895, 67.78849682660547, 14.471946008572445, 25.60008710350823, 46.07576799041101, 89.45582016906486, 34.46731622176758, 32.644767868155824, 65.79261518777068, 46.7829274144063]]
        # #
        # #
        # # 276
        # omega = [[-0.23867242270545197, -0.24963512607885785, -0.8495159563303842, 0.334571573208301, -0.14919576691032232, 0.9595608197821646, -0.5852639547534026, -0.5471656666122935, 0.021007496109239243, -0.2527090948271957, -0.8184495571242738, 0.6656316794069173, -0.8124757119245942, -0.8813543276226907, 0.39604098108724006, 0.12433338482211376, 0.851634698898561, -0.6601633314363806, 0.9639621428885572, -0.12276891567211257, 0.30956752799189924, 0.45581854756100926, -0.016480015746269894, 0.8293655911480984, 0.7015532790205834, 0.48998335955398264, -0.36624390350941427, 0.28811835997659285, 0.059219369559606605, -0.22954613784250832, 0.49877211711197433, 0.4470750439291089, 0.36054658651906224, -0.8871142116054005, -0.6519207852332063, -0.6384403361580857, 0.3188652630791573, -0.9439146069550932, -0.4738595818334199, 0.630989458733507, -0.12044047316327333, -0.23456968437785664, -0.3904946719448139, 0.242724344369029, -0.11756766478362257, 0.20913289079825303, 0.2706052209968808, -0.8910878659707913, -0.9103008276017954, -0.9965395210945227, 0.891107954739156, -0.05627011984835795, 0.598640162631618, 0.6673655475283822, -0.27654353875559146, 0.8423328363464719, 0.544144472930485, -0.01056788190269864, -0.5569438376110398, -0.8888287122695342, 0.5184076666374959, -0.8501114141060897, -0.6996939238503908, -0.7762603988928733, -0.4268673054111649, 0.13375774027049392, -0.37432473346236605, 0.691229078694712, 0.5395441084130181, -0.5937449693733752, 0.2131109287878885, 0.10917517401068522, 0.3466824946784406, -0.41426965023671136, -0.7865274686139196, 0.5431394998862937, -0.2192445631165605, 0.512238288407403, -0.08336270365828424, -0.11344445405134529, 0.13968759438730172, 0.8835526145819694, 0.8414813717069927, -0.6765354156414405, 0.0066137586382426505, 0.6453504769855953, -0.28546463848730097, -0.034912441423839535, 0.8394588394483329, -0.5635762410430276, -0.7320085971925299, -0.7885369294617974, -0.07045904939565761, -0.7619635775181484, -0.6106668332939384, 0.6236123654464225, 0.2598922136873447, 0.5371310169086123, 0.9819917326559657, 0.34843758603876696]]
        # mu = [[92.10059311702196, 6.443741268755088, 94.99999416317573, 50.655712480144246, 41.360477510827934, 54.73493648729184, 16.96830178826021, 43.34349678607248, 9.093405694330182, 5.120484745329589, 77.96400462840552, 48.72702787328448, 65.30843674006668, 14.298203151146016, 23.83284483299064, 2.286055076581439, 23.864058297354894, 15.050980684461202, 1.4240531307041393, 52.871382574863865, 83.14015154432694, 75.53547551418565, 92.2901723186651, 6.8641101295236595, 34.17222200381763, 10.284833948859495, 32.967275770400505, 61.69626530462822, 73.97521339329043, 5.228143853750167, 24.678880004044302, 39.52629303835968, 72.17786282594797, 42.43959824130822, 32.27973354969682, 19.079036147835172, 46.91501556713299, 58.96693907860965, 76.60259613799253, 97.56284615556727, 79.57399607601513, 17.268267773753344, 0.4885898686494916, 13.853635672322273, 74.34661301187055, 10.893253725375573, 12.642549852312223, 72.48417069128578, 4.282173602107775, 85.78891598388084, 93.8665230662037, 97.56661721803125, 9.083982958688496, 68.4762372606688, 17.08122012069053, 98.1287504357429, 66.64070341757606, 43.31522299857961, 87.28537495190477, 74.89110033091379, 78.66079668908034, 60.17244630401263, 31.634569243908437, 7.727893810516628, 14.119627284126047, 33.62029991110078, 41.95421820656924, 91.06977055563334, 43.52874191881406, 6.456674222161429, 56.10184190047725, 58.225055783424004, 5.397937576035494, 14.51569019494967, 22.588971184064388, 75.92499453851684, 57.131206804785684, 63.94448391927091, 60.70892646443456, 71.89125239925332, 15.719822439575438, 97.05728057456716, 12.064749752896631, 50.017038241936895, 85.78781074104765, 27.23106083697119, 45.38113451544863, 58.047434102225495, 87.58247858308043, 19.11746033753383, 0.48355327364483225, 62.233670533504984, 14.88572882755551, 59.78103592865478, 47.09237836700755, 42.82370079703999, 28.42402131903008, 40.63075705233607, 36.475770506015834, 4.301189054442212]]

        # # Both Stable outside of R
        # # 102
        #omega = [[-0.8057426858762002, -0.8920645552376902, 0.6417105660608042, 0.25582564184729617, 0.2502881309362053, -0.6933202064762034, 0.9449374754570703, -0.14090150098224363, -0.9610166389743391, -0.3455748197081394, 0.324205974890192, -0.35595583106155915, -0.3406200903612082, -0.136532065308411, 0.6651138758094712, -0.4335785826371663, -0.2463445964174038, 0.9414021292282808, -0.4394383756501048, 0.8319355230543823, 0.0857548307423659, -0.8059755426055244, 0.9992590304831026, 0.8386037594282381, 0.14219209769454388, -0.2373949421974848, -0.3110710294987944, -0.09581674596770862, 0.9363819165525027, -0.13121659794584395, 0.19750892275304932, -0.334574180308143, -0.49864841052492914, -0.6921119022734943, -0.49809776623699786, 0.6878131153664568, -0.3947584165699205, 0.9017173354345931, 0.23818181073135358, 0.3820803989530901, 0.3571484305594774, -0.9953235628864496, 0.04005299551866259, 0.591584014379053, 0.342435784672789, 0.27767862531233334, 0.9774479227263879, 0.1863335505541306, 0.4610112552922907, 0.5656138777979287, 0.36331813539682956, 0.8909890517869652, -0.49355262135227496, -0.9267121004273102, -0.9969874481672945, -0.6552160848028485, 0.6850987737554362, -0.444248245706939, 0.5669053107929882, 0.024553702820263368, -0.8000292303517191, 0.07334606949363587, -0.9337920590588875, -0.9831026927940547, -0.03193774392412063, 0.5392154867007417, -0.37013583613904544, -0.5528175050705388, 0.6690133446261637, -0.30005866559345895, -0.8651825767160655, 0.4379912841191973, -0.5901668124593253, -0.5021922973057447, 0.8011247672915649, 0.14388716610819374, 0.5641286262433978, 0.14503130088117655, -0.7423156801420778, 0.712369073168396, -0.1410616449213582, -0.13716909394115384, -0.1417257383758359, 0.9991905074089087, 0.39030784196053125, -0.16321215427132052, 0.5182454993064181, 0.03846801671906985, 0.888216706922772, -0.7435070235321819, 0.9536000513865184, 0.9934682183181669, -0.04070100408939603, 0.7659179457362781, 0.35053342418744826, 0.7686144788816254, 0.4564956144551995, -0.2420231846316605, -0.4278209761135894, -0.7066591895826455]]
        #u = [[5.084263116940146, 23.95796421141251, 38.121595920789055, 58.14073674539709, 10.654278601531642, 43.08421834459498, 34.6715259808074, 17.655798060609563, 85.47708421788899, 17.771145470278572, 36.41713601164886, 54.425895979376016, 38.11534235386061, 64.39266915228643, 32.33514246321854, 28.003174987032196, 35.68440027153832, 48.531842802649315, 61.18285985767312, 38.24483029418919, 24.40419843392614, 38.817853364248144, 12.885694126429058, 53.127747113177904, 0.5641354065358928, 31.973756519279906, 31.348655482043963, 6.733396930555447, 69.08495657634442, 15.106815581607824, 35.88858846958408, 42.65421727263667, 22.670355153164678, 54.63870820897778, 99.31750385182818, 70.15351355555563, 83.24008036640214, 36.987590534825856, 62.82994227046089, 60.57907720031399, 82.25058227910091, 48.40028098506154, 23.42375494931619, 58.35202254704927, 28.957614330289726, 36.91384866623761, 61.330300072676536, 88.90139360294089, 7.150792710146126, 81.42459270004888, 67.37670247708967, 65.70402116462301, 53.24223985007114, 53.119436612867275, 8.16411410157053, 56.37659335973081, 59.6854538452017, 72.55064791848216, 18.280976181414232, 64.69363886746736, 95.21845156888959, 68.38495324708708, 0.11719301417220107, 2.494541744633172, 8.053706278482265, 7.47811389603611, 18.403482849679488, 32.46946032534805, 95.90364997954744, 48.81367258623415, 88.81870987786237, 54.14565626485043, 56.74240462592176, 85.67184835624711, 2.4790861871840986, 32.97473901728939, 43.13486818343868, 58.01797187999183, 59.232135610755975, 70.84814412667122, 32.17219178882371, 90.80398518027496, 88.80606996724212, 57.23361714351354, 63.641650879389886, 75.09754175853807, 19.721230208083963, 74.54405724682195, 76.89774770384153, 1.2155585565759797, 61.64984938054384, 0.2782589602966401, 9.428276707566674, 48.7779879263878, 71.85427210531307, 17.41198520100513, 82.09486655855629, 3.0583631257529054, 50.623031955331264, 30.02834494786879]]
        # # 164
        #omega = [[0.39602614882386256, -0.6568439584893193, -0.4448954530547189, -0.5859143852367905, -0.49720302482639633, -0.1505111785115827, -0.4544047362437811, -0.8131404791477368, -0.8969682791948252, 0.2350692692734755, 0.6777481607917701, -0.21844183812857354, 0.8276730223933846, 0.32975441109263537, 0.5518831292411539, 0.3420286115431439, 0.16470952535595895, -0.02943697789353883, -0.9987591604825443, -0.8638684990610959, -0.1107508915079658, -0.3908342416625197, 0.5461637857581374, 0.11409909353910885, 0.5872572120030652, -0.9034839052627497, 0.15329768106941666, 0.7109812636606176, -0.44828301265480297, -0.05083827378896433, -0.42180832393042356, -0.5831138058832883, -0.7943978887790475, 0.5730416830009892, -0.8804735612330976, 0.9190613325343255, 0.3010587288976776, 0.7340023977207404, -0.41476211162971754, -0.198356050061697, 0.012746232221856113, 0.8765621052448054, 0.6866266849565867, 0.012798750948133097, -0.21574908792297487, 0.7113093086868658, -0.2849232689947918, -0.8161386789805993, -0.8689677210127089, -0.2233449164485437, 0.45410895503941306, 0.6537063089848729, -0.7520090996729387, 0.8251042963685069, -0.49521959930875314, 0.30234124516114513, -0.22168020317246317, -0.1694163789213181, 0.5253493954664945, -0.2656295018429611, 0.19347521234279452, 0.08595018545306843, -0.3683832757660155, -0.20359075593210196, -0.816479784216205, -0.6859467195648838, -0.0692665893374611, 0.2321833223594021, -0.2927136661972196, 0.9796125419016177, -0.7769699116874986, -0.447978468557229, 0.33958720826316924, 0.30211219151588686, 0.036309510849563686, 0.11481857146714547, 0.9405537513467579, 0.048310655855450246, 0.4505439536941569, -0.9204728997976475, 0.4745535396360425, -0.6031351889553711, 0.307356943713025, -0.29676118115195593, 0.9575387993076743, 0.638657938657681, -0.5894873693408444, -0.6913569308046876, 0.8831392374802653, -0.625336136470444, -0.8786940403724046, 0.7126564154743695, -0.7666025137735244, 0.5654512382180388, 0.23679871819503484, 0.49251531700156503, 0.4776529833953771, -0.39557315470246857, -0.5484667781934578, -0.9594169048667494]]
        #mu = [[0.15109367098141702, 17.304552578336185, 2.5156395988961577, 59.94887110706188, 54.467494679364606, 25.629086542925283, 45.10382636752473, 63.96247693631282, 12.580992879814145, 80.11106062422384, 62.431564585790056, 88.20355416867261, 85.73166041360012, 95.57628207110643, 10.878060407454704, 94.77028894846899, 54.852064496959066, 1.5508118251035041, 27.549611708404274, 91.17923774569073, 9.452304456449445, 90.98526961794138, 47.10220016138119, 36.58025211847901, 91.53767464330193, 6.158339704665206, 46.84105253124999, 26.480546933333958, 39.1741743600557, 21.747694500088066, 1.8988235933493969, 33.12996329414488, 23.93128543809766, 77.7268563597679, 52.797322517804965, 49.34100412177007, 95.61495384167864, 87.8851262009038, 21.561386969369913, 8.25709457562228, 23.055641564777098, 95.27865024583608, 18.157445869339583, 24.632016509989352, 92.07997077568257, 22.183748274154368, 45.973287051581636, 15.099354269038823, 99.62204553263933, 37.88425538410488, 43.962822419245576, 3.6891225474710465, 96.16964944266489, 80.66349855252868, 70.10295094666462, 74.59412552276038, 63.15811293347519, 85.13298582168454, 12.545217860875347, 47.80425613175275, 34.746335075227904, 82.47713318974557, 56.31661815496392, 87.14308319181588, 57.28831662474687, 72.86489335524816, 33.497634007477906, 17.205009971718365, 35.277654497647035, 57.107241595467805, 83.64277746512123, 61.6362383568943, 69.61728158497812, 6.494157206645246, 2.9277906985792357, 55.43145254111924, 93.5461838232726, 57.49789595596442, 95.08613878008299, 27.505382561162904, 43.58150369486478, 43.20212187482322, 4.151154474612961, 76.24090468726676, 24.991684439905892, 67.81869798093715, 11.01101227939909, 44.74684533207177, 77.64341900579382, 18.468044442033083, 15.647762558294719, 22.638062012203886, 60.6788421028788, 35.761497617839225, 50.042520302759065, 43.36731987380811, 69.55692541116404, 17.156839484596475, 35.175798354800556, 2.8420050535327146]]
        # # 193
        #omega = [[0.8813800604063355, 0.47621018906439105, 0.2695988250448962, 0.7897560205073151, 0.8244887310442128, 0.5282556815124133, 0.5957060427731795, 0.39767575568417235, 0.9873931119908217, -0.8161534720234667, -0.27207840487551804, -0.5542799280511226, 0.23651001640146552, -0.5796209094444185, -0.38765234233183676, -0.6482920305942288, -0.6949962215964329, -0.18590327391127448, 0.8346330748532349, 0.3476845467215923, -0.12285216073624872, -0.20346229787591175, 0.8096692202165565, 0.4384873096967532, -0.09717452049309006, 0.26588019680228525, 0.4948334876455349, 0.8724992210483131, 0.6388193610780841, 0.07764535051119581, -0.2937309699537596, -0.028143811476788905, 0.7550510322235529, 0.5397201953022297, -0.19658511016729974, -0.04622389429522, 0.10803423321174321, 0.8752160969536946, -0.13173413251157262, 0.25267139667117355, 0.28885816618018545, 0.8381099531685923, -0.19324357305488515, -0.4264561617728908, 0.017725174426121093, -0.19285443205111696, -0.04583296225655431, -0.3417124375678031, -0.1812837215342915, -0.8276478979141315, 0.22526831455199336, 0.7846305417592534, 0.5441459286338439, 0.3364505972984193, 0.8952536452116029, 0.19871066935394333, 0.9684396225566969, 0.700949126602822, 0.7789785279032586, -0.31748182671812986, -0.47777071540192817, -0.9673094295044302, 0.27526349747492596, -0.8339588315189781, -0.8502443601721856, 0.8581824842398809, -0.4807620894914484, 0.8318272184720705, 0.17905879108740574, 0.5436157086144431, 0.38210829821793135, 0.720588229409014, -0.6793581635270154, 0.8203315578608668, -0.6979365394850245, -0.1401748210169711, 0.7813590100815413, -0.6987492091714989, -0.8306809724579478, 0.543272657559579, -0.5944482030659477, 0.6006780434773813, -0.31321718319473724, 0.8072764121287375, -0.936252468874383, -0.46644911220139407, -0.7982407345283993, 0.3244771000206117, -0.10987599370452705, -0.6214833296691669, 0.03915851800624548, -0.9743859035192521, 0.32862201064332397, 0.3059029274281333, -0.5489115201719434, -0.4502140432256172, -0.23067552623680165, -0.03189058797506594, -0.2654446254278182, 0.6496358138867206]]
        #mu = [[64.3457146680906, 93.36298185102677, 79.97888310568345, 78.94923387117636, 66.32531546121247, 27.392884467262125, 29.557504582627757, 3.872543045139165, 21.77683340541281, 0.9973695460386844, 89.84843847366778, 2.8715437538571753, 81.91293875179774, 2.468003127197982, 58.74934021533556, 45.01092788356566, 82.53287120072468, 27.27632184299611, 50.01167121167255, 82.11183444129335, 55.66768985865307, 44.79607872528025, 38.835475328560484, 9.407402972966306, 11.87421347660248, 38.657991459636435, 13.194033526536709, 78.69751275910495, 31.27172028503792, 4.211555355802665, 75.23442236634065, 89.4775654106405, 11.421556095843744, 52.285351437130004, 71.18944291106028, 74.54661378078634, 45.59734488235644, 63.305815583894486, 45.41547833184293, 40.75717986389981, 42.439909565758704, 30.179862925637245, 80.33669513145163, 24.175059754974537, 33.56604762205504, 97.6455059378869, 6.316870387330942, 89.2308309382696, 87.93553553572401, 66.7712174046226, 0.648636615924969, 13.008045581836225, 53.67144625177289, 29.753818948261333, 32.72487913978883, 28.450917216663317, 17.204203629736924, 15.50129933962231, 31.55681757736648, 28.549410408287677, 42.11619324958629, 3.9749274343040653, 12.029598380452422, 89.10766530201263, 23.549251240054915, 84.93481224049732, 3.3729552307816157, 53.89344533139594, 82.09801638461977, 22.635847978675894, 37.337586496479055, 44.069139888576224, 2.8377280385928683, 0.7394332090028155, 47.98662522746706, 88.19700562173576, 84.50949597124591, 52.49965032186549, 96.44202914074876, 9.289572203494389, 88.26617718929208, 66.44885427029247, 51.1192164012753, 37.37094335821618, 76.17757838135248, 4.656017285835801, 12.240934380404411, 28.867624678110502, 69.29597215010855, 19.269818883162127, 36.56844980423429, 32.034873551240665, 61.723160439873034, 79.31967348060354, 70.56277094001541, 63.122059474965994, 69.36222474585463, 68.31265890871707, 29.4936960835653, 17.809338521854535]]
        #

        system_state = np.zeros(SPECIES_K+ENV_VARS)

        Eg = ENV_START[0]
        #Eg = random.uniform(0, RANGE_R)
        SURVIVAL_THRESHOLD=0.2
        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * 5 ** 2 ))
            if a_star < SURVIVAL_THRESHOLD:
                a_star = 0
            system_state[s_i] = a_star

        for _ in range(ENV_VARS):
            system_state[SPECIES_K+_] = ENV_START[_]

        ENV_VAR_ALIVE_ZERO_START = ENV_START[0]
        ENV_VAR_ALIVE_ONE_START = ENV_START[0]


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
        #plot_gaussian_trunk() - not needed
        plot_inverse_case()
        plot_temps()
        ALIVE_THRESHOLD=0.2
        plot_alphas_truncated()

        results = [[] for _ in range(SPECIES_K+ENV_VARS)]
        times_steps=[]
        system_state = np.zeros(SPECIES_K+ENV_VARS)
        Eg = ENV_START[0]

        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            if a_star < ALIVE_THRESHOLD:
                a_star = 0
            system_state[s_i] = a_star

        for _ in range(ENV_VARS):
            system_state[SPECIES_K+_] = ENV_START[_]

        for step in np.arange(TIME_START, TIME_END, TIME_STEP):
            times_steps.append(step)
            for _ in range(SPECIES_K+ENV_VARS):
                results[_].append(system_state[_])
            k1 = TIME_STEP * rates_of_change_system_state_new_niche(system_state)
            k2 = TIME_STEP * rates_of_change_system_state_new_niche(system_state + k1 * 0.5)
            k3 = TIME_STEP * rates_of_change_system_state_new_niche(system_state + k2 * 0.5)
            k4 = TIME_STEP * rates_of_change_system_state_new_niche(system_state + k3)
            system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)
        ENV_VAR_ALIVE_ONE_END = system_state[SPECIES_K+0]

        #print(results_nt[-1][-1])
        #print(results[-1][-1])

        #if((results_nt[-1][-1] > 110 or results_nt[-1][-1] < -5) and (results[-1][-1] > 110 or results[-1][-1] < -5)):
        if(1):
            #or
           # ((results_nt[-1][-1] < 100 and results_nt[-1][-1] > 0) and (results[-1][-1] > 100 or results[-1][-1] < 0))):

        #print(omega)
        #print(mu)

        # GHOST NUMBERS

        #if((results[-1][-1] > 106 or results[-1][-1] < -6)):

            #print(omega)
            #print(mu)
            #print(Eg)

            #for line in results:
            #    print(line[-2])

            print("=================================================")

            fig = plt.figure(dpi=300, figsize=(20,10))
            #fig.suptitle('Species Aliveness ' + str(sim))
            #fig.suptitle('A simulation run with 100 biotic components', fontsize=20)

            gs = fig.add_gridspec(2,2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, :])

            myList = results_nt[:-1]
            for item in myList:
                ax1.plot(times_steps,item)
            #ax1.set_title('TGAF Model', fontsize=15)
            ax1.set_xlabel('Time Steps', fontsize=12)
            ax1.set_ylabel('Abundance', fontsize=12)

            ax1.xaxis.set_minor_locator(AutoMinorLocator())
            ax1.yaxis.set_minor_locator(AutoMinorLocator())
            ax1.tick_params(which='both', width=1)
            ax1.tick_params(which='major', length=8)
            ax1.tick_params(which='minor', length=6)


            #ax1.set_ylim([0, 1])
            myList = results[:-1]
            for item in myList:
                ax2.plot(times_steps,item)
            #ax2.set_title('NW Model', fontsize=15)
            ax2.set_xlabel('Time Steps', fontsize=12)
            ax2.set_ylabel('Abundance', fontsize=12)


            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            ax2.yaxis.set_minor_locator(AutoMinorLocator())
            ax2.tick_params(which='both', width=1)
            ax2.tick_params(which='major', length=8)
            ax2.tick_params(which='minor', length=6)


            #ax2.set_ylim([0, 1])
            #ax3.set_title('System Temperature',fontsize=15)
            ax3.set_xlabel('Time Steps', fontsize=12)
            ax3.set_ylabel('Temperature', fontsize=12)
            ax3.plot(times_steps,results_nt[-1], "b", label = "TGAF Model")
            ax3.plot(times_steps, results[-1],"k", label = "IGAF Model")

            ax3.xaxis.set_minor_locator(AutoMinorLocator())
            ax3.yaxis.set_minor_locator(AutoMinorLocator())
            ax3.tick_params(which='both', width=1)
            ax3.tick_params(which='major', length=8)
            ax3.tick_params(which='minor', length=6)


            #ax3.set_ylim([0, 100])
            ax3.legend()
            plt.tight_layout()
            fig.show()
            fig.savefig(str(sim) + '.jpg')

        #number_alive_global_end = 0
            #number_alive_end = 0

            #for s_i in range(SPECIES_K):

            #    a_star = system_state[s_i]
            #    if a_star >= ALIVE_THRESHOLD:
            #        number_alive_global_end +=1

            #number_alive_end = number_alive_global_end

