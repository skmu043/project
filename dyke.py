import math, random
import matplotlib.pyplot as plt
import numpy as np


K = 100        #Number of Biotic Components
R = 100        #Essential Range (defines where Biotic Components can be present)
P = 0          #Perturbation
OE = []        #Niche
start = 0      #Time Start
end = 200      #Time End
step= 0.01      #Time Step
w = []         #Affects Parameter (ranges between -1 and 1 for each K)
u = []         #Ideal Temperature for species (between 0 and R -> the essential range)


#OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
OE = [5 for _ in range(K)]
w = [random.uniform(-1,1) for _ in range(K)]
u = [math.trunc(random.uniform(0, R)) for _ in range(K)]

N = 2           #Number of Environment Variables
E = 10          #Temperature Start value
Et = 10         #Temperature without Biotic Force

alpha = [[] for _ in range(K)] #abundance value for a species

for _ in range(K):
    alpha[_].append((math.e) ** ((-1) * (((abs(E-u[_])) ** 2) / (2*(OE[_]**2)))))


rF = []         #Biotic Force Values
rP = []         #Perturbation Values
rE = []         #Temperature with Biotic Force Values
rEt = []        #Temperature without Biotic Force Values

#Abundance values over time
rAx = [[] for x in range(K)]

time = []

biotic_force = [[] for _ in range(K)]
temperatures = []

def plot_alphas():

    for x in np.arange (0, R, step):
        temperatures.append(x)
        
    
    for y in range(K):
        for x in np.arange (0, R, step):
            biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])


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

def update(step):
    global F, P, E, Et, ra, rF, rP, rE, rEt, a, u, w
    
    fSUM = 0
    
    for _ in range(K):
        new_alpha = (math.e) ** ((-1) * (((abs((E)-u[_])) ** 2) / (2*(OE[_]**2))))

        # abundance da/dt
        alpha[_].append(alpha[_][-1] + (new_alpha - alpha[_][-1]) * step)
        rAx[_].append(alpha[_][-1])
        fSUM = fSUM + (alpha[_][-1] * w[_]) 

        # abundance directly on the graph
        #alpha[_] = al
        #rAx[_].append(alpha[_])
        #fSUM = fSUM + (alpha[_] * w[_]) # Fixed
        
    F = fSUM * 10
    P = P + (0.2 * step)
    E = E + ((P + F) * step)

    Et = Et + P

    rF.append(F)
    rP.append(P)
    rE.append(E)
    rEt.append(Et)

def plot_w():

    plt.figure(figsize=(20,5))
    plt.ylim(-1, 1)
    plt.title('Affects values for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Affects', fontsize=18)
    plt.plot(w, 'k.', label='w')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

def plot_u():
    plt.figure(figsize=(20,5))
    plt.ylim(0, R)
    plt.title('Ideal Growing Temperature for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    plt.plot(u, 'k.', label='u')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

def plot_aot():
    plt.figure(figsize=(20,10))
    plt.title('Abundance over Temperature', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for x in range(K):
        plt.plot(time,rAx[x],label = 'id %s'%x)
    plt.show()

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

def plot_efp():
    plt.figure(figsize=(30,20))
    plt.title('Simulation Values over Time', fontsize=40)
    plt.xlabel('Time', fontsize=40)
    plt.ylabel('Values', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-50, R+20)
    plt.xlim(0, end)
    plt.plot(time,rF, 'g-', label = 'biotic force')
    plt.plot(time,rP, 'b--', label = 'perturbing force(rate)')
    plt.plot(time,rE, 'r-',label = 'temperature')
    plt.plot(time,rEt, 'k.',label = 'temperature (without biota)')
    plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)
    plt.show()


for xtime in np.arange (start, end, step):
    update(step)
    time.append(xtime)
#    if(xtime == 50):
#        P += 10

#    if(xtime == 70):
#        P -= 10

plot_alphas()          #plot abundance of species over temperature
#plot_w()               #plot affects values for each species
#plot_u()               #plot ideal growing temperature for each species
plot_aot()             #plot abundance of each species over time
#plot_aot_inc_dec()     #plot species that increase temperature and decrease temperature
#plot_b_p()             #plot biotic force and P
#plot_e()               #plot temperature value over time
plot_efp()             #plot temperature, biotic force and P over time