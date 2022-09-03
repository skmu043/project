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

    #print(fXe(0.2, 5, u))

    NRange = (fXe(0.2, 5, u) - u)

    #print("NRANGE")
    #print(NRange)

    alpha = []
    temps = []

    for x in np.arange(-5, 105, 0.001):
        temps.append(x)
        alpha.append(fYaIx(x, Ni, u, NRange))
    plt.title('truncation : '+str(fYaI(fXe(0.2, 5, u), 20, u, 0)))
    plt.plot(temps,alpha)
    plt.show()


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

    for sim in range(0, 1):

        #print(sim)
        omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

        print(omega)
        print(mu)
        #omega = [[-0.7606309166232983, -0.30656459261181834, 0.05129001309671177, 0.8222774249661196, -0.542971708146222, -0.6568450067100073, -0.2496095511396157, -0.2668532256929059, 0.8881170192804011, 0.909508485820103, -0.20755931188143273, -0.49525968173430934, 0.7146094503469982, 0.858618177871677, -0.2756183587753862, 0.11248630329649356, -0.598689670721779, 0.7053711716887376, 0.28139746988367964, -0.8512075764249425, -0.7516160270573229, -0.7500594499640028, -0.9903407710753105, -0.4071967617633194, 0.9326976323171305, 0.6696594573582442, 0.5615278795832817, 0.0018042188896614242, 0.14364585275926034, -0.1475146556195548, -0.9818272182803782, -0.9606105204746622, 0.5521104337156226, -0.985524212027862, -0.4501535347912844, 0.7347080500639207, -0.48687751725120787, -0.25010603033145484, -0.7260614787652673, 0.18522993202870142, 0.2613221089767306, -0.9597428598835354, 0.49120311384107884, 0.9133916790827177, -0.5270950170393662, -0.7091124684269821, -0.06438304709440201, -0.09317320894962888, -0.9492938163564253, -0.49360237913426785, 0.698926431716375, 0.5901097772618877, -0.5179239893046841, 0.3854529991583602, 0.13240463789438417, 0.5081645866440943, 0.03433347412512022, 0.5549361572895568, -0.981220705141056, 0.7476002079344164, 0.978180008981473, -0.2037893253909986, 0.9333536806382772, 0.6266215863870039, 0.9135455527390071, -0.16301660472149537, 0.0990767640653456, -0.33110844720401644, 0.9438425879580892, 0.9274388294755305, 0.29873647062917463, 0.5817388966871686, -0.0043540073620593756, 0.49597487507407156, 0.3443365526042439, 0.1430984387907408, -0.7502893720684516, 0.36144670690924396, 0.33977413082691155, -0.5918106904314755, 0.7921032441676021, 0.04712772558509215, 0.6501512170529951, 0.3983894753380126, 0.6520404217280338, 0.018104903738153988, 0.6881594438436913, -0.05544429899984227, 0.3951868524550115, 0.1254238110056256, 0.21741232043919156, 0.894811789615579, 0.25042489666576095, -0.27455410732584595, -0.495371342773433, -0.7295573673764788, -0.3164041835355933, -0.8861918776774904, 0.5161605204331021, -0.31997066940613283]]
        #mu = [[27.47072653065391, 72.23852775955119, 16.62261751218004, 28.609485834284744, 73.01993304058304, 30.24013979458361, 57.23888786373822, 51.75784515865691, 86.23481442996228, 79.38739736743509, 26.76485401953976, 79.86784620321485, 87.19163136649745, 37.38157055724236, 60.65223516627422, 90.99013728695155, 67.94046108980127, 0.09934338840221768, 27.46472652440971, 4.442104570899874, 74.4454301381862, 43.18903024166305, 65.41425056367859, 21.772941738257035, 79.0146426820896, 73.52028427485699, 73.69158268909696, 51.65306846474306, 34.24965413161989, 44.852275128847474, 41.14957583626929, 7.773295206297581, 91.56228670579549, 62.177459233333146, 16.563610645324257, 69.36913941271399, 94.92474590221734, 10.4565993024903, 71.34167763108775, 81.18561539523029, 10.659502267218034, 59.48452232340781, 1.434446943775869, 25.935072755774545, 91.80653805141665, 77.433236458811, 16.245524578050208, 17.473934910891032, 44.96212363168001, 84.61360175038361, 45.43519383163412, 95.25995748291835, 47.592461034809794, 14.306512034645813, 68.54317884776653, 35.06874402868364, 49.9896801520091, 78.45256083336089, 80.76923451949591, 39.67799061972166, 13.037849998417927, 73.96845064664885, 35.64790567442764, 92.3768308177704, 36.548940847278125, 10.856549812549531, 92.63888936119639, 86.3373773278121, 8.091107112224648, 92.94483288818788, 62.47513902342553, 75.93483106082417, 82.00606718997655, 33.098261036676426, 22.15594427323587, 90.88174468441669, 72.72053816890872, 46.53506895978241, 74.41126995089272, 45.19184314310557, 12.482157877347388, 77.55972523656204, 91.84338411933805, 9.307824392803244, 85.37223780735305, 87.572454037445, 98.97077903829003, 15.450437287172136, 85.90015532127354, 31.68767781601015, 5.6201832599341195, 32.18538022870965, 14.539033055396033, 99.89448344574492, 84.04861338729293, 42.76573523713013, 70.01653986983355, 24.304995167303634, 96.22304536772086, 90.72238883117949]]


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
            fig.suptitle('A simulation run with 100 biotic components', fontsize=20)

            gs = fig.add_gridspec(2,2)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, :])

            myList = results_nt[:-1]
            for item in myList:
                ax1.plot(times_steps,item)
            ax1.set_title('ST Model', fontsize=15)
            ax1.set_xlabel('Time Steps', fontsize=12)
            ax1.set_ylabel('Abundance', fontsize=12)
            #ax1.set_ylim([0, 1])
            myList = results[:-1]
            for item in myList:
                ax2.plot(times_steps,item)
            ax2.set_title('ST X u=10', fontsize=15)
            ax2.set_xlabel('Time Steps', fontsize=12)
            ax2.set_ylabel('Abundance', fontsize=12)
            #ax2.set_ylim([0, 1])
            ax3.set_title('The Environment Condition',fontsize=15)
            ax3.set_xlabel('Time Steps', fontsize=12)
            ax3.set_ylabel('Temperature', fontsize=12)
            ax3.plot(times_steps,results_nt[-1], "b", label = "ST Model")
            ax3.plot(times_steps, results[-1],"k", label = "ST X Model")
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

