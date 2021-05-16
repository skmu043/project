# Python 3.9.1
import math, random
import matplotlib.pyplot as plt
import numpy as np
import time
import multiprocessing

import sys
print("Python version",sys.version)
#print("Version info.")
#print (sys.version_info)

# Fixes in this - Starting E is 0
# The essential Range starts at 0 so there can be species that will be present and some of them can bring the E down
#

K = 500       #Number of Biotic Components
R = 100        #Essential Range (defines where Biotic Components can be present)
P = 0          #Perturbation
OE = []        #Niche
start = 0      #Time Start
end = 500      #Time End
step= 0.1    #Time Step
w = []         #Affects Parameter (ranges between -1 and 1 for each K)
u = []         #Ideal Temperature for species (between 0 and R -> the essential range)

#OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
OE = [5 for _ in range(K)]
#populates affects values
w = [random.uniform(-1,1) for _ in range(K)]

while (int(len(w)) != int(len(set(w)))):
    print("Duplicate w's detected: Regenerating ...")
    w.clear()
    w = [random.uniform(-1,1) for _ in range(K)]
    
#populates optimum growing temperatures [e.g 78.213423]
u = [random.uniform(0, R) for _ in range(K)]

while(int(len(u)) != int(len(set(u)))):
    print("Duplicate u's detected: Regenerating ...")
    u.clear()
    u = [random.uniform(0, R) for _ in range(K)]

N = 2           #Number of Environment Variables
E = -100          #Temperature Start value
Et = E         #Temperature without Biotic Force

alpha = [[] for _ in range(K)] #abundance value for a species

for _ in range(K):
    alpha[_].append((math.e) ** ((-1) * (((abs((E)-u[_])) ** 2) / (2*(OE[_]**2)))))

rF = []         #Biotic Force Values
rP = []         #Perturbation Values
rE = []         #Temperature with Biotic Force Values
rEt = []        #Temperature without Biotic Force Values

#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]
#Tracks time (steps accumulation)
time = []

biotic_force = [[] for _ in range(K)]
temperatures = []

#plot abundance of species over temperature
def plot_alphas():

    for x in np.arange (0, R, step):
        temperatures.append(x)

    for y in range(K):
        for x in np.arange (0, R, step):
            biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])

    plt.figure(figsize=(25,15))
    plt.title('Biotic Force over Time', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)

    plt.show()

def biotic_alpha_parallel(_):
        global F, P, E, Et, rF, rP, rE, rEt, u, w, step

        #print(_)
        new_alpha = (math.e) ** ((-1) * (((abs((E)-u[_])) ** 2) / (2*(OE[_]**2))))
        #time scales - for each step - the next value is calculated (next abundance and next E (temperature))
        #Now with timescales in mind, simply swithcing from the current value to the newly calculated value would indicate instantaneous change
        #Instead instead of switching directly to the newly calculated value - we can approach that value via some function
        #e.g Current E=5, new E=7, instead of using E=7 we will use a function where (E=5) approaches (E=7) so the final val may be E=6

        # Keep timescales between 1 and 0 [1 = system is at the newly calculated value instantaneously whereas values closer to zero indicate slower timescales]
        # Values outside 1 and 0 will cause errors as rates would go outside model bounds
        alpha_time_scale = 1

        # abundance da/dt
        newAlpha = alpha[_][-1] + ((new_alpha - alpha[_][-1]) * step)
        alpha[_].append(alpha[_][-1] + ((newAlpha - alpha[_][-1]) * alpha_time_scale))

        rAx[_].append(alpha[_][-1])
        rAxR[_].append(alpha[_][-1] * R)
        #fSUM = fSUM + (alpha[_][-1] * w[_]) 

        #abundance directly on the graph
        #alpha[_] = al
        #rAx[_].append(alpha[_])
        #fSUM = fSUM + (alpha[_] * w[_]) # Fixed


def update(step):
    global F, P, E, Et, rF, rP, rE, rEt, u, w
    
    fSUM = 0
    
    temperature_time_scale = 0.5

    #pool = multiprocessing.Pool(processes=1)
    #pool.map(biotic_alpha_parallel, ( _ for _ in range(K)))

    #pool = multiprocessing.Pool()
    #pool.map(biotic_alpha_parallel, range(K))

    for _ in range(K):
        biotic_alpha_parallel(_)

    for _ in range(K):
        #rAx[_].append(alpha[_][-1])
        #rAxR[_].append(alpha[_][-1] * R)
        fSUM = fSUM + (alpha[_][-1] * w[_]) 

        #abundance directly on the graph
        #alpha[_] = al
        #rAx[_].append(alpha[_])
        #fSUM = fSUM + (alpha[_] * w[_]) # Fixed

    F = fSUM * 10
    P = P + (0.2 * step)
    #P = 0 
    #F = fSUM                  [Explore the linear increase for P]
    #P = P + (step/3.5)
    newE = E + ( ((P + F) * step))
    # E is old E and newE has the current value
    E = E + ((newE-E) * temperature_time_scale)

    # E not P ! This is the Temperature !
    # Incorrect one Et = Et + P
    # F becomes 0 - no biotic force as no biota
    Et = Et + ((P + 0) * step)

    rF.append(F)
    rP.append(P)
    rE.append(E)
    rEt.append(Et)

#plot affects values for each species
def plot_w():

    plt.figure(figsize=(20,5))
    plt.ylim(-1, 1)
    plt.title('Affects values for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Affects', fontsize=18)
    plt.plot(w, 'k.', label='w')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

#plot ideal growing temperature for each species
def plot_u():
    plt.figure(figsize=(20,10))
    plt.ylim(0, R)
    plt.title('Ideal Growing Temperature for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    plt.plot(u, 'k.', label='u')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

#plot abundance of each species over time where abundance is scaled up by R
def plot_aot_scaled():
    plt.figure(figsize=(20,10))
    plt.title('Abundance + Temp over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance Scaled UP + Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-50, R+20)
    plt.xlim(0, end)
    for x in range(K):
        plt.plot(time,rAxR[x],label = 'id %s'%x)

    plt.plot(time,rE, 'r.', label='E')
    plt.show()

#plot abundance of each species over time
def plot_aot():
    plt.figure(figsize=(20,10))
    plt.title('Abundance over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for x in range(K):
        plt.plot(time,rAx[x],label = 'id %s'%x)
    plt.show()

#plot species that increase temperature and decrease temperature
def plot_aot_inc_dec():
    plt.figure(figsize=(20,10))
    plt.title('Species Abundance (Blues Decrease Temperature while Reds Increase)', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for x in range(K):
        if(w[x]==0):
            plt.plot(time,rAx[x],'k-')
        if(w[x]<0):
            plt.plot(time,rAx[x],'b-')
        if(w[x]>0):
            plt.plot(time,rAx[x],'r-')
    plt.show()

#plot biotic force and P
def plot_b_p():
    plt.figure(figsize=(20,10))
    plt.title('Biotic Force and P', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Value for Biotic Force and P', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot(time,rF, 'g-', label='F')
    plt.plot(time,rP, 'b--', label='P')
    plt.show()

#plot temperature value over time
def plot_e():
    plt.figure(figsize=(20,10))
    plt.title('Environment Variable E', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot(time,rE, 'r-', label='E')
    plt.axhline(y=R)
    plt.show()

#plot temperature, biotic force and P over time
def plot_efp():
    plt.figure(figsize=(20,10))
    plt.title('Simulation Values over Time', fontsize=40)
    plt.xlabel('Time', fontsize=40)
    plt.ylabel('Values', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-80, R+20)
    plt.xlim(0, end)
    plt.plot(time,rF, 'g-', label = 'biotic force')
    plt.plot(time,rP, 'b--', label = 'perturbing force(rate)')
    plt.plot(time,rE, 'r-',label = 'temperature')
    plt.plot(time,rEt, 'k.',label = 'temperature (without biota)')
    #plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)
    plt.show()



sys.stdout.write("[%s]" % (" " * K))
sys.stdout.flush()
sys.stdout.write("\b" * (K+1))


if __name__ == '__main__':
    for xtime in np.arange (start, end, step):
        update(step)
        time.append(xtime)
        #print(xtime)

        if(xtime % 1 == 0):
            sys.stdout.write("-")
            sys.stdout.flush()
    

    sys.stdout.write("]\n")

    #plot_alphas()          #plot abundance of species over temperature
    #plot_w()               #plot affects values for each species
    #plot_u()               #plot ideal growing temperature for each species
    #plot_aot()             #plot abundance of each species over time
    plot_aot_scaled()      #plot abundance of each species over time scaled by R
    #plot_aot_inc_dec()     #plot species that increase temperature and decrease temperature
    #plot_b_p()             #plot biotic force and P
    #plot_e()               #plot temperature value over time
    plot_efp()             #plot temperature, biotic force and P over time