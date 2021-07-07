# Python 3.9.1
import math, random, sys, os, shelve, time
import matplotlib.pyplot as plt
import numpy as np


exp_name = "dyke_dual_space"
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + exp_name

# Arguments Check

if(len(sys.argv)!=13):
    print(sys.argv)
    print("Args: K, R, P, E, start, end, step, EN, OE, LP_Z, RUN_ID")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN=3, OE=5, LP_1_size, LP_2_size , RUN_ID : epoch")
    print("exit")
    sys.exit()

K = int(sys.argv[1])          #Number of Biotic Components
R = int(sys.argv[2])          #Essential Range (defines where Biotic Components can be present)
P = int(sys.argv[3])          #Perturbation
start = int(sys.argv[5])      #Time Start
end = int(sys.argv[6])        #Time End
step= float(sys.argv[7])      #Time Step
N = int(sys.argv[8])          #Number of Environment Variables
#E       = [random.uniform(0,100) for _ in range(N)]
E = [0,0,0]

F       = [0 for _ in range(N)]

ROUND = 10
#niche width
OEn = int(sys.argv[9])
OE = [OEn for _ in range(K)]

local_population_1_size = int(sys.argv[10])
local_population_2_size = int(sys.argv[11])

# Local Population 1 #############################################
local_population_1_index = []
local_population_2_index = []

for _ in random.sample(range(K-1), local_population_1_size):
    local_population_1_index.append(_)
# Local Population 1 #############################################
# Local Population 2 #############################################
for _ in range(K):
    if(_ not in local_population_1_index):
        local_population_2_index.append(_)
# Local Population 2 #############################################
#print("RandomSample")
#print(local_population_1_size)
#print(local_population_2_size)
#print("LP1: ",local_population_1_index, len(local_population_1_index))
#print("LP2: ",local_population_2_index, len(local_population_2_index))



# Local Population 2 #############################################


RUN_ID = int(sys.argv[12])


#populates affects values

# [w affects values] - for E1
# [w affects values] - for E2
# So each species has a different affects to each Environment Variable

w = [[] for _ in range(N)]
for wi in range(N):
    w[wi] = [random.uniform(-1,1) for _ in range(K)]
#print(w)
# Duplicate Check
for wi in range(N):
    while (int(len(w[wi])) != int(len(set(w[wi])))):
        print("Duplicate w's detected: Regenerating ...")
        w[wi].clear()
        w[wi] = [random.uniform(-1,1) for _ in range(K)]

#populates optimum growing temperatures [e.g 78.213423]

# [u ideal growing under E condition] - for E1
# [u ideal growing under E condition] - for E2
# So each species has a different ideal growth response to each Environment Variable

u = [[] for _ in range(N)]
for ui in range(N):
    u[ui] = [random.uniform(0, R) for _ in range(K)]

# Duplicate Check
for ui in range(N):
    while (int(len(u[ui])) != int(len(set(u[ui])))):
        print("Duplicate u's detected: Regenerating ...")
        u[ui].clear()
        u[ui] = [random.uniform(0, R) for _ in range(K)]


alpha = [[] for _ in range(K)]

# Fix applied - correctly generating alphas for global and local
for _ in range(K):
    al = []
    for ai in range(N):
        al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

    new_alpha = 0
    if(_ in local_population_1_index):
        new_alpha = (al[0] * al[1]) # If local population (product of abundances) (both local and global affect the local one)
    elif(_ in local_population_2_index):
        new_alpha = (al[0] * al[2]) # If local population (product of abundances) (both local and global affect the local one)
    else:
        new_alpha = al[0]       # Else take the first one as Eg

    alpha[_].append(new_alpha)
        #print("alpha: ",alpha)

rF = [[] for _ in range(N)]         #Biotic Force Values
rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

for _ in range(N):
    rE[_].append(E[_])                 #Input the Start Temperatures

#print("rE",rE)
#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]
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
            #print(_)
            al.append(round((math.e) ** ((-1) * (((abs((E[ei])-u[ei][_])) ** 2) / (2*(OE[_]**2)))), ROUND))
            #print("E: ", ei, " for E: ", E[ei]," with u = ", u[ei][_] , "and R = ", R,  " alpha is : ", al[-1])
            #print("E",ei," -> ",E[ei],"- a=",al[-1])
            #time scales - for each step - the next value is calculated (next abundance and next E (temperature))
            #Now with timescales in mind, simply swithcing from the current value to the newly calculated value would indicate instantaneous change
            #Instead instead of switching directly to the newly calculated value - we can approach that value via some function
            #e.g Current E=5, new E=7, instead of using E=7 we will use a function where (E=5) approaches (E=7) so the final val may be E=6
            # Keep timescales between 1 and 0 [1 = system is at the newly calculated value instantaneously whereas values closer to zero indicate slower timescales]
            # Values outside 1 and 0 will cause errors as rates would go outside model bounds


        if(_ in local_population_1_index):
            new_alpha = (al[0] * al[1]) # If local population (product of abundances) (both local and global affect the local one)
        elif(_ in local_population_2_index):
            new_alpha = (al[0] * al[2]) # If local population (product of abundances) (both local and global affect the local one)
        else:
            new_alpha = al[0]       # Else take the first one as Eg


        newAlpha = alpha[_][-1] + ((new_alpha - alpha[_][-1]) * step)
        alpha[_].append(alpha[_][-1] + ((newAlpha - alpha[_][-1]) * alpha_time_scale))

        rAx[_].append(alpha[_][-1])
        rAxR[_].append(alpha[_][-1] * R)


# The above for loop has concluded ============ new code block below ============

# Issue Detected : multiplying all Fs from both - no differentiation given to sub population

# All K affect Eg - Fg
# All sub-pop K affect only El - Fl
# F [0 is Fg and 1 is Fl same as abundance above]

    # w[0] affects EG whereas a subset of w[1] affects EL only
    for _ in range(K):
        fSUM[0] = fSUM[0] + (alpha[_][-1] * w[0][_])

    F[0] = fSUM[0] * 10
    newE = E[0] + (((0 + F[0]) * step))
    E[0] = E[0] + ((newE-E[0]) * temperature_time_scale)
    rF[0].append(F[0])
    rE[0].append(E[0])

    # ============ END EG ==============

    for _ in range(K):
        if( _ in local_population_1_index):
            fSUM[1] = fSUM[1] + (alpha[_][-1] * w[1][_])

    F[1] = fSUM[1] * 10
    newE = E[1] + (((0 + F[1]) * step))
    E[1] = E[1] + ((newE-E[1]) * temperature_time_scale)
    rF[1].append(F[1])
    rE[1].append(E[1])

    # ============ END EL1 ==============

    for _ in range(K):
        if( _ in local_population_2_index):
            fSUM[2] = fSUM[2] + (alpha[_][-1] * w[2][_])

    F[2] = fSUM[2] * 10
    newE = E[2] + (((0 + F[2]) * step))
    E[2] = E[2] + ((newE-E[2]) * temperature_time_scale)
    rF[2].append(F[2])
    rE[2].append(E[2])

    # ============ END EL2 ==============


    #for _ in range(K):
    #    for ei in range(N):
    #        fSUM[ei] = fSUM[ei] + (alpha[_][-1] * w[ei][_])

    #for ei in range(N):
    #    F[ei] = fSUM[ei] * 10
    #    newE = E[ei] + (((0 + F[ei]) * step))
    #    E[ei] = E[ei] + ((newE-E[ei]) * temperature_time_scale)
    #    rF[ei].append(F[ei])
    #    rE[ei].append(E[ei])


if __name__ == '__main__':

    E_prime             =[]
    F_prime             =[]
    alpha_prime         =[]
    rF_prime            =[]
    rE_prime            =[]
    rAx_prime           =[]
    rAxR_prime          =[]
    time_prime          =[]
    biotic_force_prime  =[]
    temperatures_prime  =[]
    simulation_run      =[]


    # First Set of calculations have occurred during initilization so appending time 0


    # sampling
    for Eg_temp in np.arange(1,100,50):
        for El_temp in np.arange(1,100,50):
            for El_temp2 in np.arange(1,100,50):
                #print("Init : ", Eg_temp, El_temp, El_temp2)
                simulation_run.append((Eg_temp,El_temp, El_temp2))
                time.append(0)
                # xtime should should start from one timestep + 0
                post_init_start = start + step
                for xtime in np.arange (post_init_start, end, step):
                    update(step)
                    time.append(xtime)

                # Going forward - after each run is done
                # Pack the data into separate data structures
                # Zero out the in-use data structures for the run
                # Re-initilize
                # Run again
                # e.g rE = [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]] >>>> rE_prime.append(rE) >>>> rE = [[] for _ in N] (the initilization bit)
                # rE_prime = [[[1,2,3,4,5,6,7][2,4,6,7,8,9,0]], [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]], [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]]]

                E_prime.append([Eg_temp, El_temp, El_temp2])
                F_prime.append(F)
                alpha_prime.append(alpha)
                rF_prime.append(rF)
                rE_prime.append(rE)
                rAx_prime.append(rAx)
                rAxR_prime.append(rAxR)
                time_prime.append(time)
                biotic_force_prime.append(biotic_force)
                temperatures_prime.append(temperatures)

                ###########################################################################################################################
                ###########################################################################################################################
                ###########################################################################################################################
                ######################################### RE INIT #########################################################################
                E = [Eg_temp,El_temp, El_temp2]
                F       = [0 for _ in range(N)]
                alpha = [[] for _ in range(K)]

                # Fix applied - correctly generating alphas for global and local
                for _ in range(K):
                    al = []
                    for ai in range(N):
                        al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

                    new_alpha = 0
                    if(_ in local_population_1_index):
                        new_alpha = (al[0] * al[1]) # If local population (product of abundances) (both local and global affect the local one)
                    elif(_ in local_population_2_index):
                        new_alpha = (al[0] * al[2]) # If local population (product of abundances) (both local and global affect the local one)
                    else:
                        new_alpha = al[0]       # Else take the first one as Eg

                    alpha[_].append(new_alpha)
                        #print("alpha: ",alpha)

                rF = [[] for _ in range(N)]         #Biotic Force Values
                rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

                for _ in range(N):
                    rE[_].append(E[_])                 #Input the Start Temperatures

                rAx = [[] for x in range(K)]
                rAxR = [[] for x in range(K)]
                time = []
                biotic_force = [[] for _ in range(K)]
                temperatures = []

                ###########################################################################################################################
                ###########################################################################################################################
                ###########################################################################################################################
                ######################################### END RE INIT #####################################################################



    plt.figure(figsize=(20,10))
    plt.title('Abundance over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)


    time = time_prime

    if(len(time[0]) > len(rAx_prime[0][0])):
        time[0].pop(0)

    for row in rAx_prime:
        for species in row:
            plt.plot(time[0],species)

    plt.savefig('aot_1.png')
    #plt.show()

    plt.figure(figsize=(20,10))
    plt.title('Abundance over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)


    abundance           = []
    abundance_local_1     = []
    abundance_local_2     = []
    abundance_not_local = []

    a_t = []
    a_l = []
    a_g = []
    a_2 = []


    for row in rAx_prime:
        for _ in range(len(time[0])):
            sum_abundance           = 0
            sum_abundance_local_1     = 0
            sum_abundance_local_2     = 0
            sum_abundance_not_local = 0

            num = 0
            for species_block in row: #(K species)
                sum_abundance += species_block[_]
                if(num in local_population_1_index):
                    sum_abundance_local_1 += species_block[num]
                elif(num in local_population_2_index):
                    sum_abundance_local_2 += species_block[num]
                else:
                    sum_abundance_not_local += species_block[num]

                num+=1
            abundance.append(sum_abundance)
            abundance_local_1.append(sum_abundance_local_1)
            abundance_local_2.append(sum_abundance_local_2)
            abundance_not_local.append(sum_abundance_not_local)

        plt.plot(time[0],abundance, linewidth=5)
        plt.plot(time[0],abundance_local_1, linewidth=5)
        plt.plot(time[0],abundance_local_2, linewidth=5)
        plt.plot(time[0],abundance_not_local, linewidth=5)

        a_t.append(abundance[-1])
        a_l.append(abundance_local_1[-1])
        a_2.append(abundance_local_2[-1])
        a_g.append(abundance_not_local[-1])

        abundance.clear()
        abundance_local_1.clear()
        abundance_local_2.clear()
        abundance_not_local.clear()
    plt.savefig('aot_2.png')
    #plt.show()

#print("total: ", a_t)
#print("total L: ", a_l)
#print("total G: ", a_g)
#print("simulation: ", simulation_run)

os.mkdir(data_directory)
# inputs used : sys.argv + date + other metadata
# temperatures, biotic_force, w, u, rAxR, time, rE, rF, rP, rE, rEt
#print(data_directory)
s = shelve.open(data_directory + "/" + exp_name + ".data")
try :
    s['sys.argv']       = sys.argv
    #s['temperatures']   = temperatures
    #s['biotic_force']   = biotic_force
    s['w']              = w
    s['u']              = u
    #s['rAx_prime']      = rAx_prime
    #s['rAxR_prime']     = rAxR_prime
    #s['time_prime']     = time_prime
    #s['rE_prime']       = rE_prime
    #s['rF_prime']       = rF_prime

    s['K']              = K
    s['R']              = R
    s['E']              = E
    #s['start']          = start
    #s['end']            = end
    #s['step']           = step
    s['N']              = N
    s['OEn']            = OEn

    s['a_t']                = a_t
    s['a_l']                = a_l
    s['a_g']                = a_g
    s['a_2']               = a_2
    s['simulation_run']     = simulation_run
    s['RUN_ID']             = RUN_ID
    s['local_population_1_size'] = local_population_1_size
    s['local_population_2_size'] = local_population_2_size
    s['local_population_1_index']  = local_population_1_index
    s['local_population_2_index']  = local_population_2_index

finally:
    s.close()

