import shelve, os
import matplotlib.pyplot as plt
import sys
import numpy as np
print(sys.version)
import math, random, time

#plot affects values for each species
def plot_w():

    plt.figure(figsize=(20,5))
    plt.ylim(-1, 1)
    plt.title('Affects values for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Affects', fontsize=18)

    for list in w:
        plt.plot(list, 'k.', label='w')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

#plot ideal growing temperature for each species
def plot_u():
    plt.figure(figsize=(20,5))
    plt.ylim(0, R)
    plt.title('Ideal Growing Temperature for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    for list in u:
        plt.plot(list, 'k.', label='u')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

#plot abundance of species over temperature
def plot_alphas():

    biotic_force = [[] for _ in range(K)]
    temperatures = []
    ROUND = 7

    s = 0.01

    for x in np.arange (-50, R+50, s):
        temperatures.append(x)

    for y in range(K):
        for x in np.arange (-50, R+50, s):
            biotic_force[y].append( (round((math.e) ** ((-1) * (((abs((x)-u[0][y])) ** 2) / (2*(5**2)))), ROUND)) * w[0][y] )

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force over Time', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(K):
        plt.plot(temperatures,biotic_force[_])
    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
    plt.show()


#plot abundance of each species over time
def plot_aot():
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

    abundance = []
    for row in rAx_prime:
        for _ in range(len(time[0])):
            sum_abundance = 0
            for species_block in row: #(K species)
                sum_abundance += species[_]
            abundance.append(sum_abundance)
        plt.plot(time[0],abundance, linewidth=5)
        abundance.clear()

    plt.show()

#plot abundance of each species over time where abundance is scaled up by R
def plot_aot_scaled():
    plt.figure(figsize=(30,30))
    plt.title('Abundance + Temp over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance Scaled UP + Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    time = time_prime

    if(len(time[0]) > len(rAx_prime[0][0])):
        time[0].pop(0)

    for row in rAxR_prime:
        # Each Row is Each Sample Run
        for species in row:
            plt.plot(time[0],species)

    for row in rE_prime:
        for env_var in row:
            plt.plot(time_prime[1],env_var, 'k-', label='E', linewidth=4)
    plt.show()

#plot temperature value over time
def plot_e():
    plt.figure(figsize=(20,10))
    plt.title('Environment Variable E', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for row in rE_prime:
        for env_var in row:
            plt.plot(time_prime[1],env_var, 'k-', label='E', linewidth=4)
    plt.axhline(y=R)
    plt.show()

#plot temperature, biotic force and P over time
def plot_efp():
    plt.figure(figsize=(30,20))
    plt.title('Simulation Values over Time', fontsize=40)
    plt.xlabel('Time', fontsize=40)
    plt.ylabel('Values', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-50, R+20)
    plt.xlim(0, end)

    time = time_prime

    if(len(time[0]) > len(rAx_prime[0][0])):
        time[0].pop(0)

    for row in rF_prime:
        for env_var in row:
            plt.plot(time[0],env_var, 'g-', label='F', linewidth=4)

    for row in rE_prime:
        for env_var in row:
            plt.plot(time_prime[1],env_var, 'r-', label='E', linewidth=4)

    plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)
    plt.show()

def plot_stable():
    #plt.figure(figsize=(30,20))
    #plt.title('Stable Points', fontsize=40)
    #plt.xlabel('E Values', fontsize=40)
    #plt.ylabel('E Values', fontsize=40)
    #plt.xticks(fontsize=20)
    #plt.yticks(fontsize=20)
    #plt.ylim(-50, R+20)
    #plt.xlim(0, end)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    #ax.view_init(0, 60)

    #ax.set_xlim3d(0,2000)
    #ax.set_ylim3d(-10,120)
    #ax.set_zlim3d(0,200)

    ax.set_xlabel('Time')
    ax.set_ylabel('EL')
    ax.set_zlabel('EG')

    for row in rE_prime:
        ax.plot(time_prime[1],row[1],row[0], label='E', linewidth=2) # rE[0] is global and goes on the y axis

    plt.show()


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):
    #print(os.listdir(data_dr))
    #print(data_dr + "/" + str(os.listdir(data_dr)[select]) + "/dyke.data.db")
    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke_space.data")
    #print(data_dr + "/" + str(data_archives[select]) + "/dyke_space.data.db")
    #print("/Users/sumeet.kumar/IdeaProjects/project/data/1624411616.2321012.dyke_space/dyke_space.data")
    #s = shelve.open("/Users/sumeet.kumar/IdeaProjects/project/data/1624411616.2321012.dyke_space/dyke_space.data")
    print(s)
    try :

        args            = s['sys.argv']
        temperatures    = s['temperatures']
        biotic_force    = s['biotic_force']
        w               = s['w']
        u               = s['u']
        time_prime      = s['time_prime']
        K               = s['K']
        R               = s['R']
        E               = s['E']
        start           = s['start']
        end             = s['end']
        step            = s['step']
        N               = s['N']
        OEn             = s['OEn']
        a_t             = s['a_t']
        a_l             = s['a_l']
        a_g             = s['a_g']
        simulation_run  = s['simulation_run']
        RUN_ID          = s['RUN_ID']

        #s['rAx_prime']      = rAx_prime
        #s['rAxR_prime']     = rAxR_prime
        #s['rE_prime']       = rE_prime
        #s['rF_prime']       = rF_prime


        while True:
            print("0 plot affects values for each species")
            print("1 plot ideal growing temperature for each species")
            print("2 plot alphas")
            print("3 plot abundance of each species over time")
            print("4 plot abundance of each species over time scaled by R")
            print("5 plot temperature value over time")
            print("6 plot temperature, biotic force")
            print("7 Stable Points")
            select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

            if select == 0:
                plot_w()               #plot affects values for each species
            elif select == 1:
                plot_u()               #plot ideal growing temperature for each species
            elif select == 2:
                plot_alphas()
            elif select == 3:
                plot_aot()             #plot abundance of each species over time
            elif select == 4:
                plot_aot_scaled()      #plot abundance of each species over time scaled by R
            elif select == 5:
                plot_e()               #plot temperature value over time
            elif select == 6:
                plot_efp()             #plot temperature, biotic force over time
            elif select == 7:
                plot_stable()        # [E1 values over time >>>>>]
                                      # [E2 values over time >>>>>]
                                      #.
                                      #.
                                      # [En values over time >>>>>]
            else:
                print("Invalid Selection")
                break
    finally:
        s.close()

else:
    print("Invalid Selection")
