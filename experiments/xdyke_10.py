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

K = 100      #Number of Biotic Components
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

#print(w)

w = [0.22290348003252802, -0.4629709973570022, -0.8305921866244728, -0.0178186215957985, -0.9669363768691881, -0.30706952478132377, -0.8761079742617301, -0.6914097743715963, 0.9489947686600657, -0.09968787770087806, 0.6654132365434, 0.5673192263992146, 0.9623053709547675, -0.9961915575180731, -0.8256564465549872, -0.7906664974173372, 0.7824980662457608, 0.8386676334388565, -0.9621026052962749, 0.05977127007991068, -0.2747672654566149, 0.8415459614629708, -0.25857959029842337, -0.4580459876488989, 0.901480494820387, 0.7562689400292626, -0.3481116519623315, 0.7113631008756001, -0.2083208377364203, -0.6365152951733213, -0.7334565356199079, 0.5972997449224304, 0.37180541308511006, -0.7586550503523757, -0.31922705389715755, 0.005823205766538386, 0.5434041409185877, 0.011155430851174541, 0.09349422839935473, -0.10846190159892832, -0.6706041334633039, -0.2712373435080153, -0.523212851682024, -0.4146555142407644, -0.7486960716430044, 0.8206772060290859, -0.42159377180171553, -0.008267104550471194, 0.9045975831671611, -0.8621166350916365, 0.04049179978248696, 0.8133262864877202, -0.7416100674566095, -0.7444340598078403, 0.07382315224109015, -0.6655606642632381, -0.9265041572201136, -0.3168412705839452, -0.543567310444161, 0.9021971815252285, -0.21690769207401184, -0.3597373220265643, -0.6167338317739857, 0.48833301685543273, -0.998134450352907, 0.2291047678679874, 0.6795725944097732, -0.5535078060463117, -0.5479218693802832, -0.7759020319758503, 0.3699969099179765, -0.19452929770266914, 0.15997709330508103, 0.22382550543000157, -0.31721663155336755, 0.8525405870041101, -0.8489099314938164, -0.22822344953026286, -0.279143977537299, -0.661111403736308, 0.5664466815267348, -0.10994570211896226, -0.5938146480075974, -0.51989552631363, 0.8002199898222355, 0.4360779413865259, 0.6622814630183638, 0.9789349384014991, -0.8339934930440012, 0.26472861028587946, 0.33444566572511647, -0.2129877234704307, -0.5624632052606644, -0.4074501004052016, -0.8838315898409441, -0.10257659539010278, 0.24480195120917725, -0.16375371920365023, 0.6502574986128222, 0.28724770888932083]
#populates optimum growing temperatures [e.g 78.213423]
u = [random.uniform(0, R) for _ in range(K)]

while(int(len(u)) != int(len(set(u)))):
    print("Duplicate u's detected: Regenerating ...")
    u.clear()
    u = [random.uniform(0, R) for _ in range(K)]

#print(u)

u = [81.54694547622039, 49.79260835911865, 76.69084216721406, 74.74185613081276, 86.37071274482399, 97.08920928660345, 65.04546515622106, 26.712647927542122, 23.68806755614685, 32.84043278269901, 98.61765125975917, 15.024820405481199, 78.75779474091067, 62.12557666746918, 94.65298920571898, 5.188318604720488, 16.744164233946314, 98.81055118345692, 25.41378250882067, 42.174106718765756, 35.69459211792299, 96.57878282579571, 44.811538140293536, 7.200627155429595, 67.79516420031784, 51.581251325510316, 18.78941781050757, 12.694540292436251, 41.60325423116772, 56.63386455847246, 65.85087069982902, 1.6425545323342616, 55.26407241819582, 2.3592759644410877, 77.1079700554468, 2.92985755961106, 52.85862986794091, 49.000412107795235, 84.8526350244889, 16.70676516036921, 76.63957828759506, 66.31003694313488, 29.42834217942274, 87.120020763632, 14.576075491389485, 40.85530244494235, 34.14766965675086, 11.146186534865897, 73.96586352428695, 75.67492337793954, 70.09157153819758, 97.71833081810058, 62.251400055610084, 93.91376421053944, 24.082459435081883, 56.37438624019566, 21.564983343628054, 12.862815935710392, 3.6534075873162153, 81.8497095183414, 66.5568333304612, 22.532826273253182, 5.360563316432332, 43.553665194832426, 21.289429188900687, 84.52475916590359, 91.54261219222016, 34.632450712033624, 90.67308757093693, 22.269761989565463, 21.61036848211595, 17.635331917612817, 22.78127897708281, 64.39452897350743, 20.949502989132352, 6.4088429323481355, 32.53204209848002, 63.73557355431853, 55.836002523919056, 91.7250626755923, 21.65964732918174, 17.74104819043182, 86.10718925405168, 48.72600740615227, 10.01939129858458, 55.91281603451155, 68.89245680369423, 92.03224594154833, 90.43429310972189, 20.492876561747142, 69.49856332174316, 24.904286316014247, 0.9104632195542361, 59.087026595151336, 14.852893198067484, 32.077622925762505, 52.757158876215385, 10.02257414035952, 83.9657699080946, 37.47411264217779]
N = 2           #Number of Environment Variables
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
    
    temperature_time_scale = 0.2

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

    F = fSUM * 15
    P = P + (0.2 * step)
    #P = 0 
    #F = fSUM                  [Explore the linear increase for P]
    #P = P + (step/3.5)
    newE = E + (((P + F) * step))
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
    #plot_aot_scaled()      #plot abundance of each species over time scaled by R
    #plot_aot_inc_dec()     #plot species that increase temperature and decrease temperature
    #plot_b_p()             #plot biotic force and P
    #plot_e()               #plot temperature value over time
    plot_efp()             #plot temperature, biotic force and P over time