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
L2 = []
L1size = []
L2size = []
L1_L2_x = [] # L1 L2 pairs - uniq

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    #print(data_archives.index(si) , si)
    #print(data_dr + "/" + str(si) + "/dyke_space.data")
    print(data_dr + "/" + str(si) + "/dyke_dual_space.data")
    s = shelve.open(data_dr + "/" + str(si) + "/dyke_dual_space.data")
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
        a_2                 = s['a_2']
        #simulation_run      = s['simulation_run']
        #RUN_ID              = s['RUN_ID']
        #local_population_   = s['local_population_']


        #print("PHI: ", args[10])
        #print(simulation_run)
        #print(RUN_ID)

        P.append(args[10])
        L1size.append(args[10]) # L1 Size
        L2size.append(args[11]) # L2 Size
        L1_L2_x.append((args[10], args[11]))
        T.append(a_t)
        G.append(a_g)
        L.append(a_l)
        L2.append(a_2)

        #s['rAx_prime']      = rAx_prime
        #s['rAxR_prime']     = rAxR_prime
        #s['rE_prime']       = rE_prime
        #s['rF_prime']       = rF_prime

    finally:
        s.close()

#print("ARGS : ",args[10],args[11])

# First get uniq items from P - essentially all the PHIs

uniq = []
for x in L1_L2_x:
    if x not in uniq:
        uniq.append(x)

uniq.sort()


Ps = []
Ts = []
Ls = []
Gs = []
L2s = []


#plt.errorbar(x, y, e, linestyle='None', marker='^')

#plt.show()

print(L1size)
print(L2size)
print(L1_L2_x)
print(uniq)

idx = 0
for phi in uniq:
    idx = 0
    for entry in L1_L2_x:
        if phi == entry:
            Ps.append(P[idx])
            Ts.append(T[idx])
            Ls.append(L[idx])
            Gs.append(G[idx])
            L2s.append(L2[idx])
        idx +=1


print(Ts)

x = np.array([])
y = np.array([])
e = np.array([])
z = []

for phi in uniq:
    print("EACH: ", phi)
    nums = []
    idx=0
    for entry in L1_L2_x:
        #print("ENTRY",entry)
        if phi == entry:
            for each_num in Ts[idx]:
                nums.append(each_num)
        idx +=1


    print(nums)
    print(phi[0], phi[1])
    x = np.append(x, str(phi[0] +"/"+ phi[1]))
    y = np.append(y, statistics.mean(nums))
    e = np.append(e, statistics.stdev(nums))
    z.append(nums)

print(x,len(x))
print(y,len(y))
print(e,len(e))

plt.figure(figsize=(20,10))
plt.title('Abundance Values at different Population Ratios', fontsize=20)
plt.xlabel('Population Ratio L1/L2', fontsize=20)
plt.ylabel('Total Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.errorbar(x, y, e, linestyle='None', marker='^', elinewidth=7, capsize=8, capthick=7)
plt.plot(x,y)
plt.plot(x, z, '.')
#plt.plot(Ps, Ts, '.' , label='A', linewidth=1)
plt.show()

plt.figure(figsize=(20,10))
plt.title('Local 1 Abundance Values at different Population Ratios', fontsize=20)
plt.xlabel('Population Ratio L1/L2', fontsize=20)
plt.ylabel('Total L1 Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(Ls, '.' , label='L1', linewidth=4)
plt.show()

plt.figure(figsize=(20,10))
plt.title('Local 2 Abundance Values at different Population Ratios', fontsize=20)
plt.xlabel('Population Ratio L1/L2', fontsize=20)
plt.ylabel('Total L2 Abundance', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.plot(L2s, '.' , label='L2', linewidth=4)
plt.show()

