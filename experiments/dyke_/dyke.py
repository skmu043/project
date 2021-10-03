# Python 3.9.1
import math, random, sys, os, shelve, time
import matplotlib.pyplot as plt
import numpy as np

print("Python version: ", sys.version, "Version info: ", sys.version_info)


exp_name = "dyke"
# the var time below conflicts with time here
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + exp_name

# Arguments Check
if(len(sys.argv)!=8):
    print("Args: K, R, P, E, start, end, step")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01")
    print("exit")
    sys.exit()

K = int(sys.argv[1])          #Number of Biotic Components
R = int(sys.argv[2])          #Essential Range (defines where Biotic Components can be present)
P = int(sys.argv[3])          #Perturbation
E = int(sys.argv[4])          #Temperature Start value
OE = []                       #Niche
start = int(sys.argv[5])      #Time Start
end = int(sys.argv[6])        #Time End
step= float(sys.argv[7])      #Time Step

Et = E          #Temperature without Biotic Force
w = []         #Affects Parameter (ranges between -1 and 1 for each K)
u = []         #Ideal Temperature for species (between 0 and R -> the essential range)

#OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
OE = [5 for _ in range(K)]
#populates affects values
w = [random.uniform(-1,1) for _ in range(K)]
#populates optimum growing temperatures
u = [math.trunc(random.uniform(0, R)) for _ in range(K)]
                #Number of Environment Variables


alpha = [[] for _ in range(K)] #abundance value for a species

for _ in range(K):
    alpha[_].append((math.e) ** ((-1) * (((abs(E-u[_])) ** 2) / (2*(OE[_]**2)))))

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
    global F, P, E, Et, rF, rP, rE, rEt, u, w
    
    fSUM = 0
    
    for _ in range(K):
        new_alpha = (math.e) ** ((-1) * (((abs((E)-u[_])) ** 2) / (2*(OE[_]**2))))

        # abundance da/dt
        alpha[_].append(alpha[_][-1] + (new_alpha - alpha[_][-1]) * step)
        rAx[_].append(alpha[_][-1])
        rAxR[_].append(alpha[_][-1] * R)
        fSUM = fSUM + (alpha[_][-1] * w[_]) 

        # abundance directly on the graph
        #alpha[_] = al
        #rAx[_].append(alpha[_])
        #fSUM = fSUM + (alpha[_] * w[_]) # Fixed

    F = fSUM * 10
    P = P + (0.2 * step)
    #F = fSUM                  [Explore the linear increase for P]
    #P = P + (step/3.5)
    E = E + ((P + F) * step)

    Et = Et + P

    rF.append(F)
    rP.append(P)
    rE.append(E)
    rEt.append(Et)

for xtime in np.arange (start, end, step):
    update(step)
    time.append(xtime)

#Create Data Dump Directory - uniq to each run

os.mkdir(data_directory)
# inputs used : sys.argv + date + other metadata
# temperatures, biotic_force, w, u, rAxR, time, rE, rF, rP, rE, rEt
print(data_directory)
s = shelve.open(data_directory + "/" + exp_name + ".data")
try :

    s['sys.argv']       = sys.argv
    s['temperatures']   = temperatures
    s['biotic_force']   = biotic_force
    s['w']              = w
    s['u']              = u
    s['rAx']              = rAx
    s['rAxR']           = rAxR
    s['time']           = time
    s['rE']             = rE
    s['rF']             = rF
    s['rP']             = rP
    s['rEt']            = rEt

    s['K'] = K
    s['R'] = R
    s['P'] = P
    s['E'] = E
    s['start'] = start
    s['end'] = end
    s['step'] = step

    s.sync()

finally:
    s.close()

