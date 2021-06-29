import shelve, os
import matplotlib.pyplot as plt
import sys
import numpy as np
print(sys.version)
import math, random, time


plt.figure(figsize=(20,10))
plt.title('Abundance Values at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Total Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

P = []
T = []
L = []
G = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)
    #print(data_dr + "/" + str(si) + "/dyke_space.data")
    s = shelve.open(data_dr + "/" + str(si) + "/dyke_space.data")
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
        #simulation_run      = s['simulation_run']
        #RUN_ID              = s['RUN_ID']
        #local_population_   = s['local_population_']


        #print("PHI: ", args[10])
        #print(simulation_run)
        #print(RUN_ID)

        P.append(args[10])
        T.append(a_t)
        L.append(a_l)
        G.append(a_g)

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
print(uniq)

Ps = []
Ts = []
Ls = []
Gs = []

for phi in uniq:
    idx = 0
    for entry in P:
        if phi == entry:
            Ps.append(P[idx])
            Ts.append(T[idx])
            Ls.append(L[idx])
            Gs.append(G[idx])
        idx +=1


plt.figure(figsize=(20,10))
plt.title('Abundance Values at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Total Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(Ps, Ts, '.' , label='A', linewidth=4)
plt.show()

plt.figure(figsize=(20,10))
plt.title('Local Abundance Values at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Total Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(Ps, Ls, '.' , label='L', linewidth=4)
plt.show()

plt.figure(figsize=(20,10))
plt.title('Global Abundance Values at different PHI Levels', fontsize=20)
plt.xlabel('PHI', fontsize=20)
plt.ylabel('Total Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(Ps, Gs, '.' , label='G', linewidth=4)
plt.show()