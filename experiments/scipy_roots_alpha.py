from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import math, random, time, sys, os, shelve
from multiprocessing import Process, Pool

#def fun(x):
#    return [x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0,0.5 * (x[1] - x[0])**3 + x[1]]
#def jac(x):
#    return np.array([[1 + 1.5 * (x[0] - x[1])**2,-1.5 * (x[0] - x[1])**2],[-1.5 * (x[1] - x[0])**2,1 + 1.5 * (x[1] - x[0])**2]])
#sol = optimize.root(fun, [0, 0], jac=jac, method='hybr')
#print(sol.x)


exp_name = "scipy_roots_alpha"
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." +exp_name

# Arguments Check
if(len(sys.argv)!=4):
    print("Args: K, N,  RUN_ID")
    print("e.g K=100, N=5, RUN_ID : epoch")
    print("exit")
    sys.exit()

K = int(sys.argv[1])          #Number of Biotic Components
N = int(sys.argv[2])          #Essential Range (defines where Biotic Components can be present)
RUN_ID = int(sys.argv[3])          #Perturbation

R = 100        #Essential Range (defines where Biotic Components can be present)
P = 0          #Perturbation
F = P
OE = []        #Niche
start = 0      #Time Start
end = 200      #Time End
step= 0.1      #Time Step
w = []         #Affects Parameter (ranges between -1 and 1 for each K)
u = []         #Ideal Temperature for species (between 0 and R -> the essential range)

#OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
OE = [N for _ in range(K)]
#populates affects values
w = [random.uniform(-1,1) for _ in range(K)]

while (int(len(w)) != int(len(set(w)))):
    print("Duplicate w's detected: Regenerating ...")
    w.clear()
    w = [random.uniform(-1,1) for _ in range(K)]

#print(w)

# #populates optimum growing temperatures [e.g 78.213423]
u = [random.uniform(0, R) for _ in range(K)]

while(int(len(u)) != int(len(set(u)))):
    print("Duplicate u's detected: Regenerating ...")
    u.clear()
    u = [random.uniform(0, R) for _ in range(K)]

#print(u)

E = -20          #Temperature Start value
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
    for x in np.arange (-30, R+30, step):
        temperatures.append(x)

    for y in range(K):
        for x in np.arange (-30, R+30, step):
            biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])

    plt.figure(figsize=(25,15))
    plt.title('Biotic Force over Time', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    for _ in range(K):
        plt.plot(temperatures,biotic_force[_])

    plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
    plt.show()

########################################################################################################################
def f1(x):
    biotic_force = []

    for y in range(K):
        biotic_force.append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[y])
    return(np.sum((np.array(biotic_force, dtype=float))))

def stable_point_return():

    x = []
    y = []

    X1 = -50
    Y1 = R + 50

    for xi in np.arange(X1, Y1, 0.1):
        x.append(xi)
        y.append(f1(xi))

    true_zeros = []
    for _ in range(R):
        sol = optimize.root(f1, [_], jac=False, method='hybr')
        if(sol.x >=0 and sol.x <= R):
            true_zeros.append(sol.x)

    zeros_uniq = []
    for xi in true_zeros:
        if(int(xi) not in zeros_uniq):
            zeros_uniq.append(int(xi))

    #print("Unique Points ...")
    zeros_uniq.sort()

    #idx=0
    #current_sign = "?"
    #if(y[idx]>0):
    #    current_sign = "+"
    #elif(y[idx]<0):
    #    current_sign = "-"

    #stable_points = []

    #for xi in np.arange(X1, Y1, 0.1):
        #print(x[idx])
    #    loopy_sign="?"
        #print(y[idx])
    #    if(y[idx]>0):
    #        loopy_sign = "+"
    #    elif(y[idx]<0):
    #        loopy_sign = "-"

    #    if(loopy_sign != current_sign):
            #print("Sign Change Detected!")
            #print("From : ", current_sign , " To : ", loopy_sign)
    #        if(current_sign == "+" and loopy_sign == "-"):
                #print(">>>>> Stable Point : ", x[idx])
    #            stable_points.append(int(x[idx]))
    #        current_sign = loopy_sign
            #print(y[idx])

    #    idx+=1

    #print(K)
    if(K == 100):
        plt.figure(figsize=(25,15))
        plt.title('Combined Biotic over Time', fontsize=40)
        plt.xlabel('Temperature', fontsize=40)
        plt.ylabel('biotic force (a * w)', fontsize=40)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.axhline(y=0, color='r', linestyle='-')
        plt.plot(x,y)
        plt.show()


    #reduced_stable_points=[]
    #for each_point in stable_points:
    #    if(each_point >=0 and each_point <= R):
    #        reduced_stable_points.append(each_point)
    new_reduced_stable_points = []

    #print(true_zeros)
    #print("----->>>>>", zeros_uniq)
    #print("Stable Points: ", stable_points)
    #print("Final Return : ", len(reduced_stable_points))

    for each_point in zeros_uniq:
        #print(each_point, ">>>" ,f1(each_point-1), int(f1(each_point)), f1(each_point+1))
        if(f1(each_point-1) > 0 and f1(each_point+1) < 0):
            new_reduced_stable_points.append(each_point)

    #print("New Stable Points", new_reduced_stable_points)

    return(len(new_reduced_stable_points))


stable_points_average = []
total_points = []

stable_point = stable_point_return()
#print("K + N + Stable : ",K , N, stable_point)

os.mkdir(data_directory)
s = shelve.open(data_directory + "/" + exp_name + ".data")
try :

    s['sys.argv']       = sys.argv
    s['w']              = w
    s['u']              = u
    s['K']              = K
    s['N']              = N
    s['stable_point']   = stable_point
    s['RUN_ID']         = RUN_ID

finally:
    s.close()
