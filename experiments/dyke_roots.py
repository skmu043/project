# Python 3.9.1
import math, random
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize


K = 100        #Number of Biotic Components
R = 100        #Essential Range (defines where Biotic Components can be present)
P = 0          #Perturbation
OE = []        #Niche
start = 0      #Time Start
end = 100      #Time End
step= 0.1     #Time Step
w = []         #Affects Parameter (ranges between -1 and 1 for each K)
u = []         #Ideal Temperature for species (between 0 and R -> the essential range)

#OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
OE = [5 for _ in range(K)]
#populates affects values
w = [random.uniform(-1,1) for _ in range(K)]
#print(w)
#w = [0.517246640841579, 0.8297022275120753, -0.05645622120783966, 0.704729947313314, 0.6524706600106398, -0.044281904417032836, -0.600116311419012, -0.34825748197661177, 0.004301944972719074, -0.9793994049648478, -0.2686037390457803, -0.08235139554740023, -0.6558644618620983, -0.1566770990905184, -0.5771731407234182, -0.7295241221399686, 0.6225845257316318, 0.9878421760516847, -0.5824535230571415, -0.31904797433881305, 0.32442502028829767, -0.42785217397706754, -0.06566644557461321, -0.014543111518315, -0.20642140135050835, 0.30129490622168476, 0.5136547975020878, -0.06359698318994877, 0.12957864225614668, 0.08969982838315183, -0.6191419975756376, -0.18139433701768404, -0.9313446173134092, -0.6678534779932797, -0.9355190572453682, -0.5053046632305578, -0.4761971897762507, 0.11611901054890694, -0.5245726178690271, -0.7427316993503972, 0.8469217827404836, 0.1724713939867175, 0.6427403306417498, -0.156361258555745, -0.16152053747207384, 0.9102633157958913, -0.5694863191368189, 0.4050481089290343, 0.40220615389920056, 0.7377159290360193, 0.6212658275919134, 0.2090209509854979, 0.0653714549794584, -0.5450818249709064, -0.6574816127804854, 0.20609405819401916, 0.7466740613009073, 0.43630090690148626, -0.7031455543038754, 0.17239753491854293, -0.0853184351845806, 0.8041335953050732, 0.2594023440882154, -0.4257742454785267, -0.12281247742922363, 0.41524848866075303, -0.19232869536434127, -0.883441515380146, -0.23993330923293255, 0.05235654811585677, 0.07583573657419662, -0.4791390943464726, 0.3551325819785631, 0.18904247390522255, 0.7229905603533391, -0.5290064707868367, -0.21685608078561258, 0.5026465536196973, -0.5267698333883888, -0.27774488331107317, -0.7582917480293556, 0.7801933046003047, -0.5980043587338015, 0.1295147087977675, 0.6746040916757079, 0.1284478651663845, -0.03442470919600171, 0.3067003918038971, 0.04281540251510618, 0.2674812018874235, 0.1671622335694758, -0.746682701603272, -0.6402295071401991, -0.879811412658227, 0.39187594004186566, 0.48170820272073, -0.9132235055348223, 0.5110207219884992, -0.8525427726678194, 0.36761149018807626]


#populates optimum growing temperatures
u = [math.trunc(random.uniform(0, R)) for _ in range(K)]
#print(u)
#u = [6, 94, 23, 9, 53, 77, 8, 4, 22, 99, 75, 99, 99, 95, 92, 81, 61, 99, 52, 12, 2, 43, 52, 13, 71, 96, 82, 63, 18, 50, 65, 11, 67, 56, 68, 86, 18, 1, 69, 33, 61, 51, 1, 98, 16, 72, 50, 85, 79, 71, 76, 64, 52, 52, 81, 95, 90, 19, 63, 7, 98, 8, 29, 88, 31, 42, 37, 20, 61, 26, 38, 13, 22, 9, 58, 62, 23, 21, 26, 1, 23, 60, 61, 17, 43, 49, 29, 89, 95, 20, 50, 87, 79, 5, 90, 98, 56, 40, 51, 44]


N = 2           #Number of Environment Variables
E = 40          #Temperature Start value
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
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]
#Tracks time (steps accumulation)
time = []

biotic_force = [[] for _ in range(K)]
temperatures = []

#plot abundance of species over temperature
def plot_alphas():

    for x in np.arange (-50, R+50, step):
        temperatures.append(x)

    for y in range(K):
        for x in np.arange (-50, R+50, step):
            biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force over Temperature', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for _ in range(K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)

    plt.show()



###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################


def f1(x):
    #return((x**3) + (2*(x**2)) - (2*x) - 5)
    #return(x**2 -1000)
    biotic_force = []
    for y in range(K):
        biotic_force.append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])

    return(np.sum((np.array(biotic_force, dtype=float))))

x = []
y = []

X1 = -50
Y1 = R + 50

for xi in np.arange(X1, Y1, 0.1):
    x.append(xi)
    y.append(f1(xi))
    print(xi," ",y[-1])

#TypeError: fsolve: there is a mismatch between the input and output shape of the 'func' argument 'f1'.Shape should be (2,) but it is (1,).

def plot_function():
    print("Plotting Sum  ... ")
    plt.figure(figsize=(20,10))
    plt.title('xy', fontsize=40)
    plt.xlabel('x', fontsize=40)
    plt.ylabel('y', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.axvline(x=0)
    plt.axhline(y=0)

    plt.plot(x,y, 'r-',label = 'roots')
    plt.show()

print("Solving ...")
true_zeros = []

for _ in range(R):
    sol = optimize.root(f1, [_], jac=False, method='hybr')
    if(sol.x >=0 and sol.x <= R):
        true_zeros.append(sol.x)

print("All points ...")
#print(true_zeros)
print("Unique Points ...")
#print(np.unique(np.array(true_zeros)))

def plot_stable_points():
    plt.figure(figsize=(20,10))
    plt.title('Roots', fontsize=40)
    plt.xlabel('temperature', fontsize=40)
    plt.ylabel('biotic force', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for stable in true_zeros:
        plt.axvline(x=stable)
    plt.axvline(x=0)
    plt.axhline(y=0)
    plt.plot(x,y, 'r-',label = 'biotic force')
    #plt.legend(loc=7, prop={'size': 30})
    plt.show()


print("Solving ...")
true_zeros = []
sign_change = ""

if(y[0] < 0):
    sign_change = "neg"
if(y[0] > 0):
    sign_change = "pos"
if(y[0] == 0):
    print("ZERO DETECTED")

print(sign_change)

for _ in range(R):
    sol = optimize.root(f1, [_], method='df-sane')
    if(sol.x >=0 and sol.x <= R):
        true_zeros.append(sol.x)

print(np.unique(np.array(true_zeros)))

###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################
###################### ROOTS ########################################################





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
    #P = P + (0.2 * step)
    P = P
    E = E + ((P + F) * step)
    # no biotic force
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
    plt.figure(figsize=(20,5))
    plt.ylim(0, R)
    plt.title('Ideal Growing Temperature for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    plt.plot(u, 'k.', label='u')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

#plot abundance of each species over time where abundance is scaled up by R
def plot_aot_scaled():
    plt.figure(figsize=(30,30))
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
    plt.figure(figsize=(30,20))
    plt.title('Simulation Values over Time', fontsize=40)
    plt.xlabel('Time', fontsize=40)
    plt.ylabel('Temperatures', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-100, R+80)
    plt.xlim(0, end)
    #plt.plot(time,rF, 'g-', label = 'biotic force')
    #plt.plot(time,rP, 'b--', label = 'perturbing force(rate)')
    plt.plot(time,rE, '.',label = 'temperature')
    #plt.plot(time,rEt, 'k.',label = 'temperature (without biota)')
    #plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)
    plt.axhline(y=0)
    #plt.axvline(x=50)
    #plt.axvline(x=55)
    plt.show()



for E_new in np.arange (-50,R+50, 10):
    E = E_new
    print("Running for E = ", E)
    for xtime in np.arange (start, end, step):
        update(step)
        time.append(xtime)

        #if(xtime==50):
        #    P += 10
        #if(xtime==55):
        #    P -= 10

plot_alphas()          #plot abundance of species over temperature

#plot_alphas()
#plot_function()
plot_stable_points()

#plot_w()               #plot affects values for each species
#plot_u()               #plot ideal growing temperature for each species
#plot_aot()             #plot abundance of each species over time
#plot_aot_scaled()      #plot abundance of each species over time scaled by R
#plot_aot_inc_dec()     #plot species that increase temperature and decrease temperature
#plot_b_p()             #plot biotic force and P
#plot_e()               #plot temperature value over time
plot_efp()             #plot temperature, biotic force and P over time


#ytime=0
#for E_new in np.arange (-50,R+50, 10):
#    E = E_new
#    print("The Line for : ", E)
#    for xtime in np.arange (start, end, step):
#        print(rE[ytime])
#        ytime+=1
#    print("---")

output = []
for x in rE:
    if math.trunc(x) not in output:
        output.append(math.trunc(x))
print(output)

print("")

for x in output:
    count = 0
    for y in rE:
        if x == math.trunc(y):
            count+=1
    if (x > 0 and x < R) and count > 100:
        print("Stable Point : ",x," > Stability : ",count)