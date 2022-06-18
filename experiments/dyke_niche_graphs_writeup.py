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
plt.rcParams["font.family"] = "Times New Roman"

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
    plt.title('Biotic Force of 100 species with alive threshold : ' + str(ALIVE_THRESHOLD), fontsize=30)
    plt.xlabel('Temperature', fontsize=20)
    plt.ylabel('Biotic Force', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(SPECIES_K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)

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
    return (u + (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))))

def fXe_negative(Ya, Ni, u):
    return (u - (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))))

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

    plt.plot(temps,alpha)
    plt.show()

    Ni = 20

    alpha = []
    temps = []

    for x in np.arange(-5, 105, 0.001):
        temps.append(x)
        alpha.append(fYaI(x, Ni, u, ST))

    plt.plot(temps,alpha)
    plt.show()

    print(fXe(0.2, 5, u))

    NRange = (fXe(0.2, 5, u) - u)

    print("NRANGE")
    print(NRange)

    alpha = []
    temps = []

    for x in np.arange(-5, 105, 0.001):
        temps.append(x)
        alpha.append(fYaIx(x, Ni, u, NRange))
    plt.title('truncation : '+str(fYaI(fXe(0.2, 5, u), 20, u, 0)))
    plt.plot(temps,alpha)
    plt.show()


    print("Truncation : " , str(fYaI(fXe(0.2, 5, u), 20, u, 0)))

    print(fXe(0.5, 20, u))
    print(fXe_negative(0.5, 20, u))

    print("Verification : ")

    print(fYaI(fXe(0.5, 20, u), 20, u, 0))
    print(fYaI(fXe_negative(0.5, 20, u), 20, u, 0))



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
    fig.suptitle('Abundance for a simulation with 20 biotic components', fontsize=30)
    #fig.set_size_inches(3, 1.5)
    for _ in range(SPECIES_K-80):
        ax1.plot(temperatures,abundance_for_st[_])
    ax1.set_title('ST Model', fontsize=20)
    ax1.set_xlabel('Temperature', fontsize=20)
    ax1.set_ylabel('Abundance', fontsize=20)
    for _ in range(SPECIES_K-80):
        ax2.plot(temperatures,abundance_for_st_x[_])
    ax2.set_title('ST X Model', fontsize=20)
    ax2.set_xlabel('Temperature', fontsize=20)
    ax2.set_ylabel('Abundance', fontsize=20)

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


    plt.figure(figsize=(20,20), dpi=300)
    #plt.title('Biotic Force for 100 species in the JI Model', fontsize=30)
    plt.xlabel('Temperature', fontsize=30)
    plt.ylabel('Biotic Force', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    for _ in range(SPECIES_K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4, label='Combined Biotic Force')
    plt.legend(prop={'size': 30})
    plt.tight_layout()
    plt.show()

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

    plt.figure(figsize=(20,20), dpi=300)
    #plt.title('Biotic Force for 100 species in the ST Model', fontsize=30)
    plt.xlabel('Temperature', fontsize=30)
    plt.ylabel('Biotic Force', fontsize=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    for _ in range(SPECIES_K):
        plt.plot(temperatures,abundance_for_st_x[_])

    plt.plot(temperatures,np.sum((np.array(abundance_for_st_x, dtype=float)), axis=0), lw=4, label='Combined Biotic Force')
    plt.legend(prop={'size': 30})
    plt.tight_layout()
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





if __name__ == '__main__':

    for sim in range(0, 1):

        print(sim)
        omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

        #omega = [[-0.45570349442317326, 0.9926566784100075, 0.888549026699764, 0.14928544806368338, -0.9996264233158934, -0.4385610485834899, -0.08038541123391796, -0.6260150923250922, -0.47235424832183703, 0.08534830008829819, 0.04489977686157176, 0.20925429233908543, 0.33254542065489434, 0.0840048379259668, -0.49418449970563194, 0.436757577051496, -0.10442494952037817, -0.4757616053069522, 0.21271766647905, -0.6493074483217449, -0.8353341996437704, -0.2056676731969611, -0.219114871126884, 0.3482756179142581, -0.7028479623473687, 0.06655864395029609, -0.2594633634515553, 0.71109331879931, 0.2747045193340545, 0.6193474305138627, -0.25427585792748153, 0.42264649771697194, -0.7670619926343019, 0.5194338631557813, 0.2955475346559173, 0.8362123688943861, 0.6573922690996239, -0.33227139325065447, -0.0417269736591952, 0.7474335021368288, 0.9721065523838919, 0.34145394392500794, -0.30158871229552675, 0.6383367973778622, -0.9658834930955591, -0.07201988587528474, 0.9859256779205499, -0.6456751377980265, 0.34596559619769796, 0.008703893399207852, -0.07077698031747182, 0.3687349388787877, -0.2128586823583034, -0.428467812969906, -0.7612417665593616, -0.753888511973025, -0.5823855135630227, 0.34720762289370155, 0.5656944867519293, 0.8403768391236615, -0.7055015387194046, 0.9167024270609578, -0.8377405059080654, 0.613368151113215, -0.20085672743369987, 0.5324877735051026, -0.5404093098424378, -0.4793030314186506, 0.020783169315794492, -0.339350185930549, 0.6054281152160448, 0.28770596557012973, -0.401032234196933, 0.38382727022568774, -0.5413820483877667, 0.758547095228723, -0.6448441963556208, 0.6560350213385717, 0.3469446743002187, 0.8625593080945051, 0.4194108163826673, -0.6339503647551701, -0.5989127516874513, 0.16883352766822046, 0.7962033825685597, 0.9763587212634572, 0.5538388523357383, -0.7229295694774867, 0.24767934527856994, -0.19051494465833296, 0.9618596596970299, -0.07639280592577236, 0.022044509842797044, -0.013989791881396485, -0.3704962409196475, 0.2672261767286863, -0.1534276130207337, -0.08114108850993196, -0.6215005516241177, 0.6759542263992802]]
        #mu = [[83.28825637081572, 19.91571692627123, 61.32525135771586, 15.672902946681667, 99.60109311987452, 24.675327758157138, 45.49766334978646, 98.01411351486603, 14.926773747336696, 26.197887508113617, 8.33039624249391, 48.66276674465474, 79.20562919962092, 31.3412946690634, 60.07521156883806, 88.51314965235781, 59.22999584861787, 10.866822137513054, 7.721663134157419, 55.212466048042465, 45.4806246853774, 85.04379254402356, 33.36264522633044, 40.3056685546258, 72.26374861765332, 89.23388249659632, 33.23888478740763, 55.56364226912627, 86.46948848547017, 60.894453999422424, 39.75123975664799, 13.596967208002864, 87.05914161816646, 49.06723745022505, 85.85137654356825, 69.49385579935682, 81.51659630137934, 66.41766460075132, 82.32530462853556, 71.40294877458498, 53.41365933231239, 72.79787798021528, 72.72213208784406, 11.622320416792142, 32.79578605848111, 77.4230920723422, 64.79175868066488, 46.11347113918881, 1.6882881330386512, 58.18262186923574, 84.53715699657215, 23.175810167449985, 33.93955774606593, 50.58356189329606, 63.62625059556708, 83.48616763999951, 92.53057854280155, 28.37362873012099, 50.22611401493394, 77.27937911974684, 9.343607616123572, 17.509062939619447, 23.933792264286257, 70.66324965971809, 9.465150222878272, 9.559033167143072, 3.9341833747434496, 11.000779992185016, 63.463144118865344, 55.13443058649155, 12.19131735262431, 40.1999622829841, 91.65020017010718, 87.98928652554477, 46.736551603767076, 67.66704553087106, 93.81008320866786, 66.1536314996456, 94.51162741284558, 69.82403754406121, 91.09060303974191, 49.9136245636434, 47.25073606399731, 60.08221965203116, 48.83208415621863, 55.86184279261995, 86.70469872121117, 4.548397985062458, 20.116473823336456, 95.11452004637383, 4.671406311441073, 48.48827517728724, 12.674309717595822, 60.89359727652629, 72.51607861597465, 1.4544709343304874, 69.00210031096566, 23.63269453291944, 96.73703677466968, 89.73427847911927]]
        #omega = [[0.7028809639906732, -0.1772091514037275, -0.5022750707740649, -0.1839283930101332, -0.573179687469179, -0.318523504945901, 0.30413821865327284, 0.04721478201458451, -0.49004990522995384, -0.9901353887100743, -0.9962422819697276, 0.5027542917212175, -0.9618118699262848, -0.7928917383683571, -0.02702868124579294, -0.9095725190225379, 0.9747470641928524, 0.6191234350155694, -0.29112650909255544, -0.5236332790394194, -0.6366900125530555, 0.35722759574639285, 0.43721823058853215, 0.9592746886641987, 0.5686812243625643, 0.9804585142478661, 0.2341045412148408, 0.46285386788226734, -0.5585416154902711, 0.01212751015140423, 0.12847071602438653, -0.6244511430766353, 0.9339738147957892, -0.8226807621014653, 0.7295867673012404, 0.17639755570437976, 0.22047086954709916, -0.15938402616898806, 0.43075611617477794, -0.2617269876194206, -0.15165831636015614, 0.07338445390080572, -0.9618530215141241, -0.395610115379051, 0.3812532804708937, -0.20232315808643153, -0.17874404398804766, 0.49894557215192314, -0.956674140516226, -0.8092590236053854, -0.6709311499649098, -0.16036194095401868, -0.5508806699696727, -0.1649714097078563, 0.7045636586948898, -0.027319879729604146, 0.4026253117955376, 0.1775493082266295, -0.5400964804660142, -0.8448081825019156, 0.5489389795162993, 0.49035250536380803, 0.8187037657352072, -0.8449921553712942, -0.12970787985505083, 0.41868472288954894, 0.3621661141351269, -0.3310855935884425, -0.33912679493185105, 0.14353678596739838, 0.036419948196754426, 0.1348853072105045, 0.831224313873234, 0.9462021869557455, 0.2300971511620058, 0.14430682431298547, -0.4958642004223879, -0.6233901784991154, 0.7278857400128496, 0.12978147103414406, 0.18889112034612054, -0.9094653448527845, -0.9866128923964586, -0.1245641678671805, 0.8351794046015812, -0.10549900467926943, -0.6288485413070777, 0.0741880417846632, -0.713221091784791, 0.9604831424854192, -0.06183853661747252, 0.4299576095696782, -0.19926841387633898, 0.10420359128516288, 0.174139571561849, 0.26640017952861106, -0.42445134222325653, 0.5032337097135968, 0.145851034601026, 0.8173483451584787]]
        #mu = [[80.74299233587666, 93.5378965976649, 10.442357422654913, 57.4125980415755, 2.6900318185959615, 65.0079767940026, 3.1695018824787424, 99.07943749361237, 95.14352033894653, 95.3905761345429, 2.827948072883979, 98.76260833009154, 54.07506651169696, 89.45656245312934, 91.88923397699847, 41.76025029034335, 56.414761402607795, 46.436050278977866, 37.496853793391615, 35.811655964306276, 26.597083761039908, 80.95536540531619, 27.55113886283106, 77.69687362612953, 75.31233676117589, 63.913588355893246, 48.96742338820174, 23.72110957687379, 1.1717537091021635, 94.15404577285385, 0.7124774563762704, 65.98862807638378, 10.252198776237897, 37.837795434272984, 9.567036513145366, 41.51553480332737, 55.23905906699075, 35.435591919211085, 87.28532523637601, 30.68815879561391, 82.29924868961325, 75.77277464255123, 58.673461463221365, 36.48404929286906, 95.3158199844202, 99.87398025953517, 25.948584989455924, 39.79944511813248, 74.49775141313818, 50.627714351469635, 2.597507513434638, 68.39765300086133, 29.689567783082282, 21.94470691531093, 37.27607963877093, 30.567639496176835, 87.78080687809955, 24.83682332506335, 81.94615990632207, 74.24416264528865, 98.47671255849522, 14.14010604180842, 45.20608600904534, 63.088617488257626, 0.20830008536402156, 23.519309768725417, 89.52238462644834, 33.34569235178243, 99.13524757044755, 40.14874380188807, 54.62058472134517, 60.691269576440654, 79.71185658177234, 99.96791311137756, 32.44832629319998, 72.71472996111216, 20.36904543494329, 32.12244788112546, 71.92141767195601, 89.7536131380911, 12.092164739612166, 31.922795916498792, 73.17438351571633, 27.76677217387672, 8.48164571815212, 91.78413911806793, 38.684544088115345, 70.66759051842355, 81.45904210444537, 38.394326951177796, 56.2409775252245, 98.4491321724531, 2.5136419711863645, 97.31579505804064, 37.20317381700997, 97.62229262211449, 46.37237645459788, 24.147839627295586, 15.552682683569419, 99.37017949594966]]


        system_state = np.zeros(SPECIES_K+ENV_VARS)

        Eg = ENV_START[0]
        #Eg = random.uniform(0, RANGE_R)

        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            system_state[s_i] = a_star

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
            k1 = TIME_STEP * rates_of_change_system_state(system_state)
            k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
            k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
            k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)
            system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)
        ENV_VAR_ALIVE_ONE_END = system_state[SPECIES_K+0]

        print(results_nt[-1][-1])
        print(results[-1][-1])

        if((results_nt[-1][-1] > 110 or results_nt[-1][-1] < -5) and (results[-1][-1] > 110 or results[-1][-1] < -5)):
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
            fig.suptitle('A simulation run with 100 biotic components', fontsize=20)

            gs = fig.add_gridspec(2,2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, :])

            myList = results_nt[:-1]
            for item in myList:
                ax1.plot(times_steps,item)
            ax1.set_title('Original Model', fontsize=15)
            ax1.set_xlabel('Time Steps', fontsize=12)
            ax1.set_ylabel('Abundance', fontsize=12)
            #ax1.set_ylim([0, 1])
            myList = results[:-1]
            for item in myList:
                ax2.plot(times_steps,item)
            ax2.set_title('Survival Threshold : ' + str(ALIVE_THRESHOLD), fontsize=15)
            ax2.set_xlabel('Time Steps', fontsize=12)
            ax2.set_ylabel('Abundance', fontsize=12)
            #ax2.set_ylim([0, 1])
            ax3.set_title('The Environment Condition',fontsize=15)
            ax3.set_xlabel('Time Steps', fontsize=12)
            ax3.set_ylabel('Temperature', fontsize=12)
            ax3.plot(times_steps,results_nt[-1], "b", label = "Original Model")
            ax3.plot(times_steps, results[-1],"k", label = "survival threshold")
            #ax3.set_ylim([0, 100])
            ax3.legend()
            fig.show()
            fig.savefig(str(sim) + '.png')

        #number_alive_global_end = 0
            #number_alive_end = 0

            #for s_i in range(SPECIES_K):

            #    a_star = system_state[s_i]
            #    if a_star >= ALIVE_THRESHOLD:
            #        number_alive_global_end +=1

            #number_alive_end = number_alive_global_end

