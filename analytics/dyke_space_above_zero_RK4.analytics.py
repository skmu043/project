import shelve, os
import matplotlib.pyplot as plt
import sys
import numpy as np
print(sys.version)
import math, random, time
import statistics


P = []
T = []
L = []
G = []

T_A     = []
T_N_A   = []
T_L_A   = []
T_L_N_A = []
T_G_A   = []
T_G_N_A = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    s = shelve.open(data_dr + "/" + str(si) + "/dyke_space_above_zero_RK4.data")
    try :
        args                = s['sys.argv']
        #temperatures        = s['temperatures']
        #biotic_force        = s['biotic_force']
        #w                   = s['w']
        #u                   = s['u']
        #time_prime          = s['time_prime']
        #K                   = s['K']
        #R                   = s['R']
        #E                   = s['E']
        #start               = s['start']
        #end                 = s['end']
        #step                = s['step']
        #N                   = s['N']
        #OEn                 = s['OEn']
        a_t                 = s['a_t']
        a_l                 = s['a_l']
        a_g                 = s['a_g']

        t_a                 = s['t_a']
        t_n_a               = s['t_n_a']
        t_l_a               = s['t_l_a']
        t_l_n_a             = s['t_l_n_a']
        t_g_a               = s['t_g_a']
        t_g_n_a             = s['t_g_n_a']
        #simulation_run      = s['simulation_run']
        RUN_ID              = s['RUN_ID']
        #local_population_   = s['local_population_']


        print("PHI: ", args[10])
        #print(simulation_run)
        #print(RUN_ID)

        P.append(args[10])
        T.append(a_t)
        L.append(a_l)
        G.append(a_g)

        T_A.append(t_a)
        T_N_A.append(t_n_a)
        T_L_A.append(t_l_a)
        T_L_N_A.append(t_l_n_a)
        T_G_A.append(t_g_a)
        T_G_N_A.append(t_g_n_a)

        #s['rAx_prime']      = rAx_prime
        #s['rAxR_prime']     = rAxR_prime
        #s['rE_prime']       = rE_prime
        #s['rF_prime']       = rF_prime

    finally:
        s.close()


# First get uniq items from P - essentially all the PHIs

uniq = []
for x in P:
    if x not in uniq:
        uniq.append(x)

uniq.sort()

################
Ps = [] # Phi
Ts = [] # Total Abundance
Ls = [] # Local Abundance
Gs = [] # Global Abundance
###############
Zs = [] # Above Zero Abundance
###############
T_As     = []
T_N_As   = []
T_L_As   = []
T_L_N_As = []
T_G_As   = []
T_G_N_As = []
##############

#plt.errorbar(x, y, e, linestyle='None', marker='^')

#plt.show()

for phi in uniq:
    idx = 0
    for entry in P:
        if phi == entry:
            Ps.append(P[idx])
            Ts.append(T[idx])
            Ls.append(L[idx])
            Gs.append(G[idx])
            T_As.append(T_A[idx])
            T_N_As.append(T_N_A[idx])
            T_L_As.append(T_L_A[idx])
            T_L_N_As.append(T_L_N_A[idx])
            T_G_As.append(T_G_A[idx])
            T_G_N_As.append(T_G_N_A[idx])
        idx +=1

# This Set Does Total Abundance Sum, Mean and Standard Deviation (Error Bar Plot)
x = np.array([])
y = np.array([])
e = np.array([])

# This Set Does Abundance above Zero Totals Sum, Mean and Standard Deviation (Error Bar Plot)
zx = np.array([])
zy = np.array([])
ze = np.array([])

# This Set Does Abundance above Zero for Locals Sum, Mean and Standard Deviation (Error Bar Plot)
zlx = np.array([])
zly = np.array([])
zle = np.array([])

# This Set Does Abundance above Zero for Globals Sum, Mean and Standard Deviation (Error Bar Plot)
zgx = np.array([])
zgy = np.array([])
zge = np.array([])


for phi in uniq:
    nums    = []
    znums   = []
    zlnums  = []
    zgnums  = []
    idx=0
    for entry in Ps:
        if phi == entry:
            # extracts exact abundance value for Totals
            #print(Ts[idx])
            for each_num in Ts[idx]:
                nums.append(each_num)
            # extracts non zero abundance values for Totals
            for each_num in T_As[idx]:
                znums.append(each_num)
            # extracts non zero abundance values for Locals
            for each_num in T_L_As[idx]:
                zlnums.append(each_num)
            # extracts non zero abundance values for Globals
            for each_num in T_G_As[idx]:
                zgnums.append(each_num)

        idx +=1

    x = np.append(x, phi)
    y = np.append(y, statistics.mean(nums))
    e = np.append(e, statistics.stdev(nums))

    zx = np.append(zx, phi)
    zy = np.append(zy, statistics.mean(znums))
    ze = np.append(ze, statistics.stdev(znums))

    zlx = np.append(zlx, phi)
    zly = np.append(zly, statistics.mean(zlnums))
    zle = np.append(zle, statistics.stdev(zlnums))

    zgx = np.append(zgx, phi)
    zgy = np.append(zgy, statistics.mean(zgnums))
    zge = np.append(zge, statistics.stdev(zgnums))

plt.figure(figsize=(20,10))
plt.title('Abundance Values at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Total Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.errorbar(x, y, e, linestyle='None', marker='^', elinewidth=7, capsize=8, capthick=7)
plt.plot(x,y)
plt.plot(Ps, Ts, '.' , label='A', linewidth=1)
plt.show()
plt.savefig(str(str(RUN_ID) + "_1.png"))

plt.figure(figsize=(20,10))
plt.title('Number of Alive Species at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Above Zero Number Alive', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.errorbar(zx, zy, ze, linestyle='None', marker='^', elinewidth=7, capsize=8, capthick=7)
plt.plot(zx,zy)
plt.plot(Ps, T_As, '.' , label='A', linewidth=1)
plt.show()
plt.savefig(str(str(RUN_ID) + "_2.png"))

plt.figure(figsize=(20,10))
plt.title('Number of Local Species Alive at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Above Zero Number Local Alive', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.errorbar(zlx, zly, zle, linestyle='None', marker='^', elinewidth=7, capsize=8, capthick=7)
plt.plot(zlx,zly)
plt.plot(Ps, T_L_As, '.' , label='A', linewidth=1)
plt.show()
plt.savefig(str(str(RUN_ID) + "_3.png"))

plt.figure(figsize=(20,10))
plt.title('Number of Global Species Alive at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Above Zero Number Global Alive', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.errorbar(zgx, zgy, zge, linestyle='None', marker='^', elinewidth=7, capsize=8, capthick=7)
plt.plot(zgx,zgy)
plt.plot(Ps, T_G_As, '.' , label='A', linewidth=1)
plt.show()
plt.savefig(str(str(RUN_ID) + "_4.png"))