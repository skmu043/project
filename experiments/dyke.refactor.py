import sys
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#mpl.use('macosx')

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

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']

finally:
    s.close()


SAMPLE_STEP = 50

K = biotic_components_K
R = essential_range_R
P = external_perturbation_rate_P
start = time_start
end = time_end
step = time_step
N = environment_components_N
E = [0, 0]
F = biotic_force_F

ROUND = truncated_gaussian_ROUND

OEn = niche_width
OE = [OEn for _ in range(K)]

local_population_ = local_population_size
local_population_index = []
local_population_size = int(local_population_/100 * K)

uniq_k = []
for x in range(local_population_size):
    one = random.randint(0,K-1)
    while one in uniq_k:
        one = random.randint(0,K-1)
    uniq_k.append(one)
    local_population_index.append(one)


local_population_index.sort()

w = affects_w
u = optimum_condition_u

#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]

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

        #print("alpha: ",alpha)

rF = [[] for _ in range(N)]         #Biotic Force Values
rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

for _ in range(N):
    rE[_].append(E[_])                 #Input the Start Temperatures

#print("rE",rE)

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

        # al = [0.5,0.5] (each species has its abundance calculated for E1 and E2 based on individual u values)
        new_alpha = 0
        if( _ in local_population_index):
            new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
        else:
            new_alpha = al[0]       # Else take the first one as Eg

        # For each species : If Species

        k1 = new_alpha
        k2 = k1 + (k1 * step/2)
        k3 = k1 + (k2 * step/2)
        k4 = k1 + (k3 * step)
        yt = alpha[_][-1] + (((k1 + (2*k2) + (2*k3) + k4)/6) - alpha[_][-1]) * step

        alpha[_].append(yt)

        #rAx[_].append(alpha[_][-1])
        #rAxR[_].append(alpha[_][-1] * R)
        rAx[_].append(yt)
        rAxR[_].append(yt * R)

    # ALPHA is the CHANGE to the previous Alpha !


    # The above for loop has concluded ============ new code block below ============

    # Issue Detected : multiplying all Fs from both - no differentiation given to sub population

    # All K affect Eg - Fg
    # All sub-pop K affect only El - Fl
    # F [0 is Fg and 1 is Fl same as abundance above]

    # w[0] affects EG whereas a subset of w[1] affects EL only
    for _ in range(K):
        fSUM[0] = fSUM[0] + (alpha[_][-1] * w[0][_])

    F[0] = fSUM[0] #* 10                #FINALLY 10 has been removed !
    newE = E[0] + (((0 + F[0]) * 1))    # [* 1] changes to [* step] -> if [* 10] is brought back above
    E[0] = E[0] + ((newE-E[0]) * temperature_time_scale)
    rF[0].append(F[0])
    rE[0].append(E[0])

    # ============ END EG ==============

    for _ in range(K):
        if( _ in local_population_index):
            fSUM[1] = fSUM[1] + (alpha[_][-1] * w[1][_])

    F[1] = fSUM[1] #* 10                #FINALLY 10 has been removed !
    newE = E[1] + (((0 + F[1]) * 1))    # [* 1] changes to [* step] -> if [* 10] is brought back above
    E[1] = E[1] + ((newE-E[1]) * temperature_time_scale)
    rF[1].append(F[1])
    rE[1].append(E[1])

    # ============ END EL ==============


if __name__ == '__main__':


    E_prime             =[]
    F_prime             =[]
    alpha_prime         =[]
    rF_prime            =[]
    rE_prime            =[]
    rAx_prime           =[]
    rAxR_prime          =[]
    rAxS_prime          =[]
    time_prime          =[]
    biotic_force_prime  =[]
    temperatures_prime  =[]
    simulation_run      =[]


    # First Set of calculations have occurred during initilization so appending time 0


    # sampling
    for Eg_temp in np.arange(1,100,SAMPLE_STEP):
        for El_temp in np.arange(1,100,SAMPLE_STEP):
            print("Init : ", Eg_temp, El_temp)
            simulation_run.append((Eg_temp,El_temp))
            time.append(0)
            # xtime should should start from one timestep + 0
            post_init_start = start + step
            for xtime in np.arange (post_init_start, end, step):
                update(step)
                time.append(xtime)
            #rAx.insert(0,0)
            #rAxR.insert(0,0)

            # Going forward - after each run is done
            # Pack the data into separate data structures
            # Zero out the in-use data structures for the run
            # Re-initilize
            # Run again
            # e.g rE = [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]] >>>> rE_prime.append(rE) >>>> rE = [[] for _ in N] (the initilization bit)
            # rE_prime = [[[1,2,3,4,5,6,7][2,4,6,7,8,9,0]], [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]], [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]]]

            E_prime.append([Eg_temp, El_temp])
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
            E = [Eg_temp,El_temp]
            F       = [0 for _ in range(N)]
            alpha = [[] for _ in range(K)]

            # Fix applied - correctly generating alphas for global and local
            for _ in range(K):
                al = []
                for ai in range(N):
                    al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

                    new_alpha = 0
                    if( _ in local_population_index):
                        new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
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


def plot_alphas():

    for x in np.arange (-50, R+50, step):
        temperatures.append(x)

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force over Temperature', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    for each_env_var in range(N):
        biotic_force = [[] for _ in range(K)]
        for y in range(K):
            for x in np.arange (-50, R+50, step):
                biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[each_env_var][y])) ** 2) / (2*(OE[y]**2)))) * w[each_env_var][y])


        for _ in range(K):
            plt.plot(temperatures,biotic_force[_])

        plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
    plt.show()



plot_alphas()


def stable_points_space():

    stable_locations = []

    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Regions', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    for row in rE_prime:

        if((int(row[1][-1]),int(row[0][-1])) not in stable_locations):
            stable_locations.append((int(row[1][-1]),int(row[0][-1])))

        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            # float color must be between 0 and 1
            # trajectories outside 0 R
            plt.plot(row[1][0], row[0][0], marker='.', markersize = "10", color=(float(0), float(0), float(1))) # Plots Start but not the Ends
            plt.plot(row[1],row[0], label='E', linewidth=1, color=(float(0), float(0), float(1)))
        else:
            plt.plot(row[1],row[0], label='E', linewidth=1, color=(float(c_r/100), float(c_g/100), float(0.5)))

            plt.plot(row[1][0], row[0][0], marker='.', markersize = "10" , color=(float(c_r/100), float(c_g/100), float(0.5)))
            plt.plot(row[1][-1], row[0][-1], marker='*', markersize = "10" , color=(float(c_r/100), float(c_g/100), float(0.5)))

    #plt.savefig("tra_reg_rgb" + str(RUN_ID) + "-" + str(random.randint(100, 999)) +".png")
    plt.show()

stable_points_space()

def stable_points_space_3d_rotate():

    fig = plt.figure(figsize=(50,50))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Total Abundance")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Abundance', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])

    index_A = 0

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_abundance = []
        row_abundance.append(0)
        decompose = rAx_prime[index_A]
        #print(decompose)
        run_length = len(decompose[0])
        #print(run_length)


        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_abundance.append(current_sum)

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_abundance[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], 0, marker='.', s=10, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_abundance[-1], marker='*', s=15, color=(float(c_r/100), float(c_g/100), float(0.5)))

    index_A +=1

    #plt.savefig("3d_abundance_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

stable_points_space_3d_rotate()