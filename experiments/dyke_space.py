# Python 3.9.1
import math, random, sys, os, shelve, time
import matplotlib.pyplot as plt
import numpy as np


exp_name = "dyke_space"
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + exp_name

# Arguments Check
if(len(sys.argv)!=10):
    print("Args: K, R, P, E, start, end, step, EN, OE")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN=2, OE=5")
    print("exit")
    sys.exit()

K = int(sys.argv[1])          #Number of Biotic Components
R = int(sys.argv[2])          #Essential Range (defines where Biotic Components can be present)
P = int(sys.argv[3])          #Perturbation
start = int(sys.argv[5])      #Time Start
end = int(sys.argv[6])        #Time End
step= float(sys.argv[7])      #Time Step
N = int(sys.argv[8])          #Number of Environment Variables
E       = [random.uniform(0,100) for _ in range(N)]
print("E's : ", E)
F       = [0 for _ in range(N)]

#niche width
OEn = int(sys.argv[9])
OE = [OEn for _ in range(K)]

#populates affects values

# [w affects values] - for E1
# [w affects values] - for E2
# So each species has a different affects to each Environment Variable

w = [[] for _ in range(N)]
for wi in range(N):
    w[wi] = [random.uniform(-1,1) for _ in range(K)]

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

# Future Fix - this calculates abundance for only the first Environment Variable

for _ in range(K):
    al = []
    for ai in range(N):
        al.append((math.e) **
                  ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))))
    alpha[_].append(np.prod(al))


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

#plot abundance of species over temperature
def plot_alphas():

    if N == 1:
        for x in np.arange (0, R, step):
            temperatures.append(x)

        for y in range(K):
            for x in np.arange (0, R, step):
                biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[0][y])) ** 2) / (2*(OE[y]**2)))) * w[y])

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

    else:
        print("N is greater than 1 - cannot visualize higher dimensions .... yet :)")


def update(step):
    global F, P, E, Et, rF, rP, rE, rEt, u, w, N

    fSUM = [0 for _ in range(N)]
    alpha_time_scale = 0.7
    temperature_time_scale = 0.2

    for _ in range(K):
        al = []
        for ei in range(N):
            #print(_)
            al.append((math.e) ** ((-1) * (((abs((E[ei])-u[ei][_])) ** 2) / (2*(OE[_]**2)))))
            #time scales - for each step - the next value is calculated (next abundance and next E (temperature))
            #Now with timescales in mind, simply swithcing from the current value to the newly calculated value would indicate instantaneous change
            #Instead instead of switching directly to the newly calculated value - we can approach that value via some function
            #e.g Current E=5, new E=7, instead of using E=7 we will use a function where (E=5) approaches (E=7) so the final val may be E=6

            # Keep timescales between 1 and 0 [1 = system is at the newly calculated value instantaneously whereas values closer to zero indicate slower timescales]
            # Values outside 1 and 0 will cause errors as rates would go outside model bounds
        new_alpha = 0
        # abundance da/dt


        if( _ in local_population_index):
            new_alpha = np.prod(al)
        else:
            new_alpha = al[0]

        newAlpha = alpha[_][-1] + ((new_alpha - alpha[_][-1]) * step)
        alpha[_].append(alpha[_][-1] + ((newAlpha - alpha[_][-1]) * alpha_time_scale))

        rAx[_].append(alpha[_][-1])
        rAxR[_].append(alpha[_][-1] * R)

    for _ in range(K):
        for ei in range(N):
            fSUM[ei] = fSUM[ei] + (alpha[_][-1] * w[ei][_])


    for ei in range(N):
        F[ei] = fSUM[ei] * 10

        # P is 0
        newE = E[ei] + (((0 + F[ei]) * step))
        # E is old E and newE has the current value
        E[ei] = E[ei] + ((newE-E[ei]) * temperature_time_scale)

        # E not P ! This is the Temperature !
        # Incorrect one Et = Et + P
        # F becomes 0 - no biotic force as no biota
        #Et = Et + ((P[ei] + 0) * step)

        rF[ei].append(F[ei])
        rE[ei].append(E[ei])


if __name__ == '__main__':
    #print("rE",rE)
    for population_size_percent in np.arange(0 , 100 , 10):
        print(population_size_percent,"%")
        local_population_size = int(population_size_percent/100 * K)
        local_population_index = []
        for x in range(local_population_size):
            local_population_index.append(random.randint(0,K-1))

        for xtime in np.arange (start, end, step):
            update(step)
            time.append(xtime)

        E       = [random.uniform(0,100) for _ in range(N)]
        for _ in range(N):
            rE[_].append(E[_])                 #Input the Start Temperatures

        print("E's : ", E)
    #print("rE",rE)
        #if(xtime % 1 == 0):
            #    sys.stdout.write("-")
            #    sys.stdout.flush()
        #print("")

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

    s['K'] = K
    s['R'] = R
    s['E'] = E
    s['start'] = start
    s['end'] = end
    s['step'] = step
    s['N'] = N
    s['OEn'] = OEn

finally:
    s.close()

