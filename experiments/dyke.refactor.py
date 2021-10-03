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

print(w)
print(u)

#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]

rNumberAlive = [[] for x in range(K)]

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

        if(new_alpha > 0):
            rNumberAlive[_].append(1)
        else:
            rNumberAlive[_].append(0)

print("alpha: ",alpha)
print("Number Alive : ", rNumberAlive)

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
    E[0] = E[0] + ((newE-E[0]) * temperature_time_scale)
    rF[0].append(F[0])
    rE[0].append(E[0])

    # ============ END EG ==============

    for _ in range(K):
        if( _ in local_population_index):
            fSUM[1] = fSUM[1] + (alpha[_][-1] * w[1][_])

    F[1] = fSUM[1]
    newE = E[1] + (((0 + F[1]) * 1))
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
    rNumberAlive_prime  =[]


    # First Set of calculations have occurred during initilization so appending time 0


    # sampling
    for Eg_temp in np.arange(1,100,S_STEP):
        for El_temp in np.arange(1,100,S_STEP):
            print("Init : ", Eg_temp, El_temp)
            simulation_run.append((Eg_temp,El_temp))
            time.append(0)
            # xtime should should start from one timestep + 0
            post_init_start = start + step
            for xtime in np.arange (post_init_start, end, step):
                update(step)
                time.append(xtime)

            #print("time: ", time)
            print("Ltime", len(time))
            E_prime.append([Eg_temp, El_temp])
            F_prime.append(F)
            for _ in range(K):
                del alpha[_][-1]
            alpha_prime.append(alpha)
            print("abundance : ",len(alpha[1]))
            rF_prime.append(rF)
            print("rF: ",len(rF[0]))
            print("rF: ", len(rF[1]))
            rE_prime.append(rE)
            print("rE: ",len(rE[0]))
            print("rE: ",len(rE[1]))
            for _ in range(K):
                 del rAx[_][-1]
            rAx_prime.append(rAx)
            print("rAx: ", len(rAx[0]))
            for _ in range(K):
                del rAxR[_][-1]
            rAxR_prime.append(rAxR)
            print("rAxR: ", len(rAxR[0]))
            for _ in range(K):
                del rNumberAlive[_][-1]
            rNumberAlive_prime.append(rNumberAlive)
            print("rNumberAlive: ",len(rNumberAlive[0]))
            time_prime.append(time)
            biotic_force_prime.append(biotic_force)
            temperatures_prime.append(temperatures)

            #fig = plt.figure(figsize=(50,50))
            #ax = fig.add_subplot(111, projection='3d', adjustable='box')
            #ax.set_title(label = "EL/EG with Number Alive Single")
            #ax.set_xlabel('X - EL', fontsize=10)
            #ax.set_ylabel('Y - EG', fontsize=10)
            #ax.set_zlabel('Z - Total Alive', fontsize=10)

            #c_r = int(rE[1][-1])
            #c_g = int(rE[0][-1])

            #print(len(rE[0]))
            #print(len(rE[1]))
            #print(len(rNumberAlive))
            #print(len(rNumberAlive[0]))
            #sum_alives = []
            #for x in range (len(rNumberAlive[0])):
            #    current_sum = 0
            #    for _ in range(K):
            #        current_sum += rNumberAlive[_][x]
            #    sum_alives.append(current_sum)

            #print(len(sum_alives))
            #if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            #    ax.scatter(rE[1][0], rE[0][0], sum_alives[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #    plt.plot(rE[1],rE[0], sum_alives , color=(float(0), float(0), float(1)))
            #else:
            #    plt.plot(rE[1],rE[0],sum_alives, color=(float(c_r/100), float(c_g/100), float(0.5)))
            #    ax.scatter(rE[1][0], rE[0][0], 0, marker='.', s=10, color=(float(c_r/100), float(c_g/100), float(0.5)))
            #    ax.scatter(rE[1][-1], rE[0][-1], sum_alives[-1], marker='*', s=15, color=(float(c_r/100), float(c_g/100), float(0.5)))

            #plt.show()

            ###########################################################################################################################
            ###########################################################################################################################
            ###########################################################################################################################
            ######################################### RE INIT #########################################################################
            E = [Eg_temp,El_temp]
            F       = [0 for _ in range(N)]
            alpha = [[] for _ in range(K)]
            rAx = [[] for x in range(K)]
            rAxR = [[] for x in range(K)]
            rNumberAlive = [[] for x in range(K)]

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
                    rAx[_].append(new_alpha)
                    rAxR[_].append(new_alpha * R)
                    if(new_alpha > 0):
                        rNumberAlive[_].append(1)
                    else:
                        rNumberAlive[_].append(0)

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

    super_biotic_force = []

    for each_env_var in range(N):
        biotic_force = [[] for _ in range(K)]
        for y in range(K):
            for x in np.arange (-50, R+50, step):
                biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[each_env_var][y])) ** 2) / (2*(OE[y]**2)))) * w[each_env_var][y])


        for _ in range(K):
            plt.plot(temperatures,biotic_force[_])


        plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
        super_biotic_force.append(np.sum((np.array(biotic_force, dtype=float)), axis=0))

    sum = []
    #for _ in range(time_end/time_step):
    #    sum.append(super_biotic_force[0][_] + super_biotic_force[1][_])
    #plt.plot(temperatures, sum, lw=10)
    #plt.show()
#plot_alphas()


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

mpl.use('macosx') #for the 3D magic

def stable_points_space_3d_rotate():

    fig = plt.figure(figsize=(300,300))
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

        print(len(row[1]))
        print(len(row[0]))
        print(len(rAx_prime[index_A]))


        row_abundance = []
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
            ax.scatter(row[1][0], row[0][0], row_abundance[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_abundance[-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))

        index_A +=1

    #plt.savefig("3d_abundance_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

stable_points_space_3d_rotate()


def stable_points_space_3d_rotate_number_alive():

    fig = plt.figure(figsize=(300,300))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Number Alive")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Alive', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])

    index_A = 0

    #print(rNumberAlive_prime)

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_alive = []

        decompose = rNumberAlive_prime[index_A]
        run_length = len(decompose[0])

        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_alive.append(current_sum)

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_alive[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_alive , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_alive, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], row_alive[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_alive[-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))

        index_A +=1

    #plt.savefig("3d_alive_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

stable_points_space_3d_rotate_number_alive()