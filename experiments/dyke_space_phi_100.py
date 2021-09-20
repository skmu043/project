# Python 3.9.1
import math, random, sys, os, shelve, time
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import matplotlib as mpl
mpl.use('macosx')
#mpl.use('TkAgg')

print(mpl.rcsetup.interactive_bk)
print(mpl.rcsetup.non_interactive_bk)
print(mpl.rcsetup.all_backends)

exp_name = "dyke_space"
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
# str(random.randint(100, 999)) added above to remove the condition where directory is not created as it exists (scheduler running at the same time step)

# Arguments Check
if(len(sys.argv)!=12):
    print("Args: K, R, P, E, start, end, step, EN, OE, LP_Z, RUN_ID")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN=2, OE=5, LP_Z = (10 - 100), RUN_ID : epoch")
    print("exit")
    sys.exit()
SAMPLE_STEP=50
K = int(sys.argv[1])          #Number of Biotic Components
R = int(sys.argv[2])          #Essential Range (defines where Biotic Components can be present)
P = int(sys.argv[3])          #Perturbation
start = int(sys.argv[5])      #Time Start
end = int(sys.argv[6])        #Time End
step= float(sys.argv[7])      #Time Step
N = int(sys.argv[8])          #Number of Environment Variables
#E       = [random.uniform(0,100) for _ in range(N)]
E = [0,0]

F       = [0 for _ in range(N)]

ROUND = 10
#niche width
OEn = int(sys.argv[9])
OE = [OEn for _ in range(K)]

local_population_ = int(sys.argv[10])
# Local Population #############################################
local_population_index = []
local_population_size =  int(local_population_/100 * K)

uniq_k = []
for x in range(local_population_size):
    one = random.randint(0,K-1)
    while(one in uniq_k):
        one = random.randint(0,K-1)
    uniq_k.append(one)
    local_population_index.append(one)


local_population_index.sort()
#print(local_population_index)

# print("Local Population Index : ", local_population_index)
# Local Population #############################################


RUN_ID = int(sys.argv[11])


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
print(w)
#w = [[-0.692497649427541, -0.7245792619554263, 0.6169230012989777, 0.6786537400709256, -0.09455620035006662, 0.5559148564241729, 0.7164217803143125, 0.4013236672189364, -0.14178255882057567, -0.9747605756370852, 0.40966825325408096, 0.6934147639849109, 0.7736733051705806, 0.3035667154951851, -0.8995273572267259, 0.5290568413158803, 0.04216859015763852, 0.0019720658302697647, -0.07002152256326344, 0.7479809598056926, 0.12295553059278896, -0.928523998501104, -0.6071817634092247, 0.48724452134440366, -0.4310893062834702, 0.7333010699160389, 0.6908312214165944, -0.30264189986423484, 0.7485891976167327, -0.9917139317318426, -0.14291326261348702, 0.031980211967836514, -0.4298048738074871, -0.11304958467032433, 0.06507279495479912, -0.8715641173583555, 0.3071049717049752, 0.6086396162742214, -0.954454028032286, -0.40355646678116197, 0.7400150002323871, 0.7546272069723237, 0.6019326363296129, -0.5995177667590821, 0.9940316267282205, 0.47463592725111203, -0.19877977490117127, 0.5726095530191737, -0.41531091923192176, 0.30029258255616265, -0.7501693441739508, -0.26865703162467813, 0.8200857678453268, 0.636087454530005, 0.9795063194602467, 0.28085592605109166, -0.08696085829032096, -0.768010239666207, 0.5557526899053955, -0.01287340733885678, 0.09649602106256316, -0.15082161859632803, 0.024268168692763714, -0.679451623652221, -0.35067464991651187, -0.7990864313551154, 0.5456788687462439, 0.6260151479683602, -0.09444349510814254, -0.049301972966594665, 0.0026825817400681906, 0.12263609575594492, 0.28708371007572, -0.5584154958869951, 0.8425889533817141, -0.26079495611396686, 0.19011583551382571, 0.5642130023185585, 0.5748210630975581, -0.9496337656382696, 0.4038653525826106, -0.5099754843165305, -0.03517570856607066, -0.5783815710029716, 0.6507452452884566, 0.25557107514775446, -0.15659017276277498, -0.8309439361641704, -0.10164017201285636, 0.02851935216668089, 0.48700003250195567, -0.5737078336111778, 0.4359604731603488, -0.9110841935010268, -0.3519630939113454, 0.514575868756953, -0.11874305677870578, -0.8688246862809201, -0.5303478151365801, -0.8528610794566827], [0.8224181036394356, -0.7864454511615318, -0.46843234047510296, -0.6792101903898924, 0.4749475040114568, 0.5261696271416121, 0.44597624693875226, -0.8801941493669709, -0.037003968866697434, 0.9189674906210215, -0.6089328135365504, -0.7604566904996459, -0.7369609016958154, -0.2978907505776216, -0.7349269961038627, -0.8674846980537052, 0.11807814037687803, -0.3030470648764225, 0.041719933020136546, 0.8080051336539245, -0.19169365395499582, 0.267319343549117, -0.33926101460610636, 0.11741401641106708, -0.898941378637466, -0.38303326683901817, -0.5370410682288589, 0.8216717947305139, -0.6016284955603282, 0.7228510461785325, -0.45139928346944447, 0.03891765201731934, -0.8018618789391183, 0.6721364473479032, -0.5495125279360618, -0.671462383098206, -0.6794545687703044, 0.46369592605301, 0.09375580839856301, -0.24741101131543597, -0.36237450074516686, 0.8548760801546949, -0.966223339203282, -0.2103064740285132, -0.6828767226870827, 0.653478091832312, 0.9398474632011173, 0.28136277563362855, 0.9760759118234839, 0.7368706742399058, 0.6405959131638688, 0.055893878893287896, -0.9240876990597895, 0.9062616390319498, 0.05335669878383653, 0.47987977145223315, 0.69113546637083, 0.053038873936661846, 0.4140686122598176, -0.2126856561876822, 0.23816840312778842, -0.12405656473772719, 0.6776272962531131, -0.06659043319474489, -0.07951029663028275, -0.09535091207085156, 0.27447270207355956, 0.6397986396578095, -0.967758200880035, -0.7683516208605305, 0.43832800533028116, -0.48749617825005176, -0.26873814516656935, 0.8502310500269423, 0.3381178221453969, -0.1126549976518807, 0.6589678422232914, -0.7638435474751089, 0.9737194730532435, 0.49484455307647246, 0.8662709314221255, 0.25987530452248087, 0.04045645823469668, 0.19210046263672886, -0.3023372330211811, -0.13262631016245185, 0.18656944649346063, -0.29349615289473396, -0.6722589749502292, -0.16676101047593295, -0.29536564487770334, 0.03322182906102533, 0.49886328266883084, -0.9176925527582696, 0.3909280344697206, 0.5198334859558658, -0.7285279842942614, 0.7008918323958304, 0.6496073245192002, 0.1453499163562728]]

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

print(u)
#u = [[65.5066228949152, 11.835087065028006, 23.648589287488917, 69.00789422527114, 51.15657020624074, 52.442972507085884, 49.82692703747995, 69.17272193331284, 67.81083081436272, 13.843638039273653, 13.038133378484906, 50.4757906454788, 51.47181254308336, 88.64914171645434, 12.195931565573591, 35.435138902020434, 52.064724203482896, 75.36915904602184, 49.18737759370373, 12.318366570398132, 24.018541559811958, 12.786424507243254, 62.863855415657625, 93.34735373175953, 57.64964145618252, 8.678329540370978, 37.03305239315996, 98.0996741403074, 72.28701025828987, 2.5750440324715185, 26.28563172960647, 23.963541765674044, 57.31133996948691, 63.762219171817414, 23.091280526009093, 8.088820603154566, 28.88122346588753, 55.85379976078709, 26.944873994365924, 52.15199541247364, 83.5700044498731, 31.25895529110172, 48.778920150837166, 94.10995636916238, 97.53089760669926, 8.587613284329121, 83.87745499324927, 18.81146009450826, 94.69804642874239, 26.44130774834146, 7.6746126266344135, 3.84069411077006, 98.70607183819797, 69.31695826472851, 55.17340352566779, 83.80515135428126, 65.17814928124217, 87.75443575593332, 55.89257806035227, 24.69999878319109, 41.71684547424959, 9.112949365976197, 81.61953528410162, 58.453359687045406, 20.64981321654107, 40.50034615375851, 43.51164596162478, 68.86638174432575, 46.888949142237514, 56.45644738171225, 34.614304536351895, 45.55601831202206, 21.829776259107135, 24.73347701283598, 36.16158486564169, 94.25533121070262, 31.52676914080227, 78.30826134003482, 56.52173678832057, 71.30644085249202, 41.86775684308392, 21.805860020502234, 56.99712790191671, 70.9983702121563, 53.26713056321376, 32.63932144634537, 97.83013952519877, 13.67912266043625, 69.97189048662045, 70.665281966703, 15.536828807169167, 53.1202404930679, 80.86392672763778, 8.660634144467993, 71.75110962002375, 62.46156336786798, 51.570723652647196, 61.10369056443429, 58.751327336547945, 76.84239020924319], [7.894449958559013, 82.81383333656042, 71.37344623728701, 48.68882119112764, 73.09093540379216, 31.097278583836196, 59.03685193456174, 67.97067265888288, 77.15908230848481, 51.60238880373922, 47.80538745241708, 79.97158482281796, 11.050635398971009, 29.193846979378158, 17.026654324120592, 32.04065305749375, 49.725586791386256, 43.99027927483633, 59.42748969694779, 27.878235303166974, 11.75078705370205, 68.78938770067792, 26.394763274907064, 52.03002483728532, 19.384481711099667, 20.8502876268714, 74.51108513996626, 56.444648630720614, 73.59412807952937, 41.39485387023083, 97.17087162493073, 32.74738976541867, 66.46321872129852, 68.40783858081267, 32.59815031260111, 96.46205544906492, 25.894693084075836, 80.89740695261482, 4.370552998947397, 8.85904722714519, 23.46071630359383, 2.0379732906211534, 43.268772482978704, 4.4232862598378775, 43.51077088193658, 30.04354063669197, 84.70492322660982, 59.98775989020948, 35.567340817727136, 63.83573914002064, 1.8092215584734506, 15.979480550014703, 80.1864737769951, 92.43685770086331, 67.31886171779252, 82.2817898191472, 84.64132702082183, 27.266895710642558, 93.47023616098944, 5.0326107245271645, 12.517820901028752, 10.405244036387506, 45.93912113077688, 69.70452199151205, 83.56092622172694, 22.51968813930658, 70.00713156062562, 2.059136490601021, 89.06568344573998, 10.892313182813208, 91.9441163829679, 53.221397452591034, 85.27019723598987, 39.02148582193424, 64.30871026045352, 66.44658345607031, 60.971583576786024, 32.824252015425095, 7.667314710266149, 45.826524417585425, 11.038768654544807, 55.817057092002806, 42.01411839609669, 33.560746524496054, 2.967832490629252, 82.78654696367043, 55.94827318544079, 92.38003985674075, 7.056471602681691, 83.15881242137549, 34.90683390856917, 21.277425207483713, 63.93015315794889, 83.10671138848242, 63.71084738686533, 21.131516877819223, 22.586195467256985, 70.14337915289308, 77.12511391935635, 98.80671700510345]]

#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]

alpha = [[] for _ in range(K)]

# Fix applied - correctly generating alphas for global and local
for _ in range(K):
    al = []
    for ai in range(N):
        al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

        new_alpha = 0
        if( _ in local_population_index):
            new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
        else:
            new_alpha = al[0]       # Else take the first one as Eg

        alpha[_].append(new_alpha)


        #print("alpha: ",alpha)

rF = [[] for _ in range(N)]         #Biotic Force Values
rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

for _ in range(N):
    rE[_].append(E[_])                 #Input the Start Temperatures

#print("rE",rE)

#Tracks time (steps accumulation)
time = []

biotic_force = [[] for _ in range(K)]
temperatures = []


def update(step):
    global F, P, E, Et, rF, rP, rE, rEt, u, w, N

    fSUM = [0 for _ in range(N)]
    alpha_time_scale = 0.7
    temperature_time_scale = 0.2

    for _ in range(K):
        al = []
        for ei in range(N):
            #print(_)
            al.append(round((math.e) ** ((-1) * (((abs((E[ei])-u[ei][_])) ** 2) / (2*(OE[_]**2)))), ROUND))
            #print("E: ", ei, " for E: ", E[ei]," with u = ", u[ei][_] , "and R = ", R,  " alpha is : ", al[-1])
            #print("E",ei," -> ",E[ei],"- a=",al[-1])
            #time scales - for each step - the next value is calculated (next abundance and next E (temperature))
            #Now with timescales in mind, simply swithcing from the current value to the newly calculated value would indicate instantaneous change
            #Instead instead of switching directly to the newly calculated value - we can approach that value via some function
            #e.g Current E=5, new E=7, instead of using E=7 we will use a function where (E=5) approaches (E=7) so the final val may be E=6
            # Keep timescales between 1 and 0 [1 = system is at the newly calculated value instantaneously whereas values closer to zero indicate slower timescales]
            # Values outside 1 and 0 will cause errors as rates would go outside model bounds

        # al = [0.5,0.5] (each species has its abundance calculated for E1 and E2 based on individual u values)
        new_alpha = 0
        if( _ in local_population_index):
            new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
        else:
            new_alpha = al[0]       # Else take the first one as Eg

        # For each species : If Species

        k1 = new_alpha
        k2 = k1 + (k1 * step/2)
        k3 = k1 + (k2 * step/2)
        k4 = k1 + (k3 * step)
        yt = alpha[_][-1] + (((k1 + (2*k2) + (2*k3) + k4)/6) - alpha[_][-1]) * step

        alpha[_].append(yt)

        #rAx[_].append(alpha[_][-1])
        #rAxR[_].append(alpha[_][-1] * R)
        rAx[_].append(yt)
        rAxR[_].append(yt * R)

# ALPHA is the CHANGE to the previous Alpha !


# The above for loop has concluded ============ new code block below ============

# Issue Detected : multiplying all Fs from both - no differentiation given to sub population

# All K affect Eg - Fg
# All sub-pop K affect only El - Fl
# F [0 is Fg and 1 is Fl same as abundance above]

    # w[0] affects EG whereas a subset of w[1] affects EL only
    for _ in range(K):
        fSUM[0] = fSUM[0] + (alpha[_][-1] * w[0][_])

    F[0] = fSUM[0] #* 10                #FINALLY 10 has been removed !
    newE = E[0] + (((0 + F[0]) * 1))    # [* 1] changes to [* step] -> if [* 10] is brought back above
    E[0] = E[0] + ((newE-E[0]) * temperature_time_scale)
    rF[0].append(F[0])
    rE[0].append(E[0])

    # ============ END EG ==============

    for _ in range(K):
        if( _ in local_population_index):
            fSUM[1] = fSUM[1] + (alpha[_][-1] * w[1][_])

    F[1] = fSUM[1] #* 10                #FINALLY 10 has been removed !
    newE = E[1] + (((0 + F[1]) * 1))    # [* 1] changes to [* step] -> if [* 10] is brought back above
    E[1] = E[1] + ((newE-E[1]) * temperature_time_scale)
    rF[1].append(F[1])
    rE[1].append(E[1])

    # ============ END EL ==============


if __name__ == '__main__':

    E_prime             =[]
    F_prime             =[]
    alpha_prime         =[]
    rF_prime            =[]
    rE_prime            =[]
    rAx_prime           =[]
    rAxR_prime          =[]
    rAxS_prime          =[]
    time_prime          =[]
    biotic_force_prime  =[]
    temperatures_prime  =[]
    simulation_run      =[]


    # First Set of calculations have occurred during initilization so appending time 0


    # sampling
    for Eg_temp in np.arange(1,100,SAMPLE_STEP):
        for El_temp in np.arange(1,100,SAMPLE_STEP):
            print("Init : ", Eg_temp, El_temp)
            simulation_run.append((Eg_temp,El_temp))
            time.append(0)
            # xtime should should start from one timestep + 0
            post_init_start = start + step
            for xtime in np.arange (post_init_start, end, step):
                update(step)
                time.append(xtime)
            #rAx.insert(0,0)
            #rAxR.insert(0,0)

            # Going forward - after each run is done
            # Pack the data into separate data structures
            # Zero out the in-use data structures for the run
            # Re-initilize
            # Run again
            # e.g rE = [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]] >>>> rE_prime.append(rE) >>>> rE = [[] for _ in N] (the initilization bit)
            # rE_prime = [[[1,2,3,4,5,6,7][2,4,6,7,8,9,0]], [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]], [[1,2,3,4,5,6,7][2,4,6,7,8,9,0]]]

            E_prime.append([Eg_temp, El_temp])
            F_prime.append(F)
            alpha_prime.append(alpha)
            rF_prime.append(rF)
            rE_prime.append(rE)
            rAx_prime.append(rAx)
            rAxR_prime.append(rAxR)
            time_prime.append(time)
            biotic_force_prime.append(biotic_force)
            temperatures_prime.append(temperatures)

            ###########################################################################################################################
            ###########################################################################################################################
            ###########################################################################################################################
            ######################################### RE INIT #########################################################################
            E = [Eg_temp,El_temp]
            F       = [0 for _ in range(N)]
            alpha = [[] for _ in range(K)]

            # Fix applied - correctly generating alphas for global and local
            for _ in range(K):
                al = []
                for ai in range(N):
                    al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

                    new_alpha = 0
                    if( _ in local_population_index):
                        new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
                    else:
                        new_alpha = al[0]       # Else take the first one as Eg

                    alpha[_].append(new_alpha)
                    #print("alpha: ",alpha)

            rF = [[] for _ in range(N)]         #Biotic Force Values
            rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

            for _ in range(N):
                rE[_].append(E[_])                 #Input the Start Temperatures

            rAx = [[] for x in range(K)]
            rAxR = [[] for x in range(K)]
            time = []
            biotic_force = [[] for _ in range(K)]
            temperatures = []

            ###########################################################################################################################
            ###########################################################################################################################
            ###########################################################################################################################
            ######################################### END RE INIT #####################################################################


def aot():
    plt.figure(figsize=(20,10))
    plt.title('Abundance over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    time = time_prime

    if(len(time[0]) > len(rAx_prime[0][0])):
        time[0].pop(0)

    for row in rAx_prime:
        for species in row:
            plt.plot(time[0],species)

    plt.show()

def aot_l_g():

    plt.figure(figsize=(20,10))
    plt.title('Local vs Global Abundance over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)


    abundance           = []
    abundance_local     = []
    abundance_not_local = []

    a_t = []
    a_l = []
    a_g = []


    for row in rAx_prime:
        for _ in range(len(time[0])):
            sum_abundance           = 0
            sum_abundance_local     = 0
            sum_abundance_not_local = 0

            num = 0
            for species_block in row: #(K species)
                sum_abundance += species_block[_]
                if(num in local_population_index):
                    sum_abundance_local += species_block[num]
                else:
                    sum_abundance_not_local += species_block[num]

                num+=1
            abundance.append(sum_abundance)
            abundance_local.append(sum_abundance_local)
            abundance_not_local.append(sum_abundance_not_local)

        plt.plot(time[0],abundance, linewidth=5)
        plt.plot(time[0],abundance_local, linewidth=5)
        plt.plot(time[0],abundance_not_local, linewidth=5)

        a_t.append(abundance[-1])
        a_l.append(abundance_local[-1])
        a_g.append(abundance_not_local[-1])

        abundance.clear()
        abundance_local.clear()
        abundance_not_local.clear()

    plt.show()


    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

def stable_points_space_kmeans():

    ### NEEDS COMPLETING IF NEEDED - not fully implemented !!! ###

    plt.figure(figsize=(30,20))
    plt.title('Trajectories', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    #fig = plt.figure(figsize=(10,10))
    #ax = fig.add_subplot(111, projection='3d')
    #ax.view_init(0, 0)
    #ax.set_axis_off()
    #ax.set_xlim3d(0,2000)
    #ax.set_ylim3d(-10,120)
    #ax.set_zlim3d(0,200)

    #ax.set_xlabel('Time')
    #ax.set_ylabel('EL')
    #ax.set_zlabel('EG')

    stable_locations = []

    for row in rE_prime:
        #ax.plot(time_prime[1],row[1],row[0], label='E', linewidth=2) # rE[0] is global and goes on the y axis
        plt.plot(row[1],row[0], label='E', linewidth=2) # rE[0] is global and goes on the y axis

        #print("EGs : ", row[0][0], "EGe : ", row[0][-1])
        #print("ELs : ", row[1][0], "ELe : ", row[1][-1])
        # x, y
        if((int(row[1][-1]),int(row[0][-1])) not in stable_locations):
            stable_locations.append((int(row[1][-1]),int(row[0][-1])))

    plt.show()

    plt.figure(figsize=(30,20))
    plt.title('Regions', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)


    print("Stable Locations: ", stable_locations)
    valid_stable_locations = []
    for stable_point_xy in stable_locations:
        if(stable_point_xy[0] >= 0 and stable_point_xy[1] >= 0 and stable_point_xy[0] <= K and stable_point_xy[1] <= K):
            valid_stable_locations.append((stable_point_xy[0], stable_point_xy[1]))

    print("Valid Stable Locations: ", valid_stable_locations)

    ############### KMeans Clusters #################################
    # Gets Number of Clusters
    sil = []
    kmax = int(len(valid_stable_locations))

    for k in range(2, kmax):
        kmeans = KMeans(n_clusters = k).fit(valid_stable_locations)
        labels = kmeans.labels_
        sil.append(silhouette_score(valid_stable_locations, labels, metric = 'euclidean'))

    print("Sil Vals : ", sil)

    max_value = max(sil)
    max_index = sil.index(max_value)
    print("Number of Clusters Needed for Kmeans : ", max_index+1)

    # Finds the points from the number of clusters defined above

    kmeans = KMeans(n_clusters=int(max_index+1), init='k-means++', max_iter=300, n_init=10, random_state=0)
    pred_y = kmeans.fit_predict(valid_stable_locations)
    plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], s=300, c='red')
    plt.show()


    cluster_i = 0

    reduced_stable_points = []

    for _ in range(len(kmeans.cluster_centers_[:, 0])):
        reduced_stable_points.append((int(kmeans.cluster_centers_[:, 0][_]), int(kmeans.cluster_centers_[:, 1][_])))
        cluster_i+=1

    region_plots_x = [[] for _ in range(len(valid_stable_locations))]
    region_plots_y = [[] for _ in range(len(valid_stable_locations))]
    region_plots_x_ = [[] for _ in range(len(valid_stable_locations))]
    region_plots_y_ = [[] for _ in range(len(valid_stable_locations))]


    print("Reduced Stable Point : ", reduced_stable_points)

    stable_index = 0

    reduced_stack = [([],[]) for _ in len(reduced_stable_points)]

    for item in reduced_stable_points:
        reduced_stack.append((item,[]))

    current_shortest_distance = 1000
    current_shortest_distance_xy = (0,0)



    # Notes : We get valid_stable_locations > these are anything in the 0 - 100 0 - 100 grid : anything venturing off to the outside are not tracked
    # now with higher sample points, and increased times, clusters form i.e they stabilize at points around a region not converging to a single point
    # that's ok so with that -> we are doing KMeans to get clusters -> and their centers so we can graph lines converging to that cluster
    # to converge to that point (same color type things)
    # reduced_stable_points are the coordinates for those points

    # Pending work : map valid_stable_locations to reduced_stable_points
    # Produce :
    # [
    # (reduced_stable_point), [(valid_stable_location), (valid_stable_location), (valid_stable_location)],
    # (reduced_stable_point), [(valid_stable_location)],
    # (reduced_stable_point), [(valid_stable_location), (valid_stable_location), (valid_stable_location)],
    # (reduced_stable_point), [(valid_stable_location), (valid_stable_location)]
    # ]
    # and then graph, for each valid_stable_location, graph to its corresponding reduced_stable_point for very high numbers



    for each_stable_point in valid_stable_locations:
        for each_reduced_stable_point in reduced_stable_points:
            #/ (x2 - x1) + (y2 - y1)
            dist = (math.sqrt( (each_stable_point[1] - each_reduced_stable_point[0]) + (each_stable_point[0] - each_reduced_stable_point[1]) ))
            if ( dist <  current_shortest_distance ):
                current_shortest_distance = dist
                current_shortest_distance_xy = (each_reduced_stable_point[0], each_reduced_stable_point[1])
        current_shortest_distance = 1000
        reduced_stack_i = 0
        for item in reduced_stack:
            if(item[0] == current_shortest_distance_xy):
                reduced_stack[reduced_stack_i].append(current_shortest_distance_xy)
            reduced_stack_i += 1


    for row in rE_prime:
        if(int(row[1][-1]) == each_stable_point[0] and int(row[0][-1]) == each_stable_point[1]):
                region_plots_x[stable_index].append(row[1])
                region_plots_y[stable_index].append(row[0])
                region_plots_x_.append(row[1][-1])
                region_plots_y_.append(row[0][-1])




    stable_index +=1
    print("Each Stable Point Done")
    colors = ['r','g','c','m','b','k','y']
    color_i = 0

    print("Valid Stable Locations : ", valid_stable_locations)
    print("Valid X : ", region_plots_x_)
    print("Valid Y : ", region_plots_y_)



    for _ in range(len(reduced_stable_points)):
        if(color_i == 7):
            color_i = 0
        plt.plot(region_plots_x[_], region_plots_y[_], '.', color = colors[color_i], markersize="1", zorder=0)
        color_i +=1
    print("Plotting Regions Done")
    for _ in range(len(reduced_stable_points)):

        #plt.scatter(valid_stable_locations[_][0],valid_stable_locations[_][1], s=180, facecolors='none', edgecolors='k', zorder=10)
        plt.scatter(reduced_stable_points[_][0],reduced_stable_points[_][1], s=180, facecolors='none', edgecolors='k', zorder=10)

    plt.show()

    #plt.plot(row[1],row[0], label='E', linewidth=2) # rE[0] is global and goes on the y axis


    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################



def stable_points_space():

    #plt.figure(figsize=(30,30))
    #plt.title('Trajectories', fontsize=40)
    #plt.xlabel('EL', fontsize=40)
    #plt.ylabel('EG', fontsize=40)
    #plt.xticks(fontsize=20)
    #plt.yticks(fontsize=20)
    #plt.ylim(-20, R+20)
    #plt.xlim(-20, R+20)

    stable_locations = []

    #for row in rE_prime:
    #    plt.plot(row[1],row[0], label='E', linewidth=2) # rE[0] is global and goes on the y axis

    #    if((int(row[1][-1]),int(row[0][-1])) not in stable_locations):
    #        stable_locations.append((int(row[1][-1]),int(row[0][-1])))
    #plt.savefig("tra_reg_vanilla" + str(RUN_ID) + str(random.randint(100, 999)) +".png")
    #plt.show()


    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Regions', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    ### HOW THIS IS DONE ###
    # EG EL GRID -> 100 x 100
    # RGB value for each refined (0.1 or 0.01) point on the grid
    # The trajectory gets the color of which ever mini RGB grid it arrives at
    # Idea by Matthew

    # 3
    # 2
    # 1
    # 0 1 2 3 .... 100

    #for row in rE_prime:
    #    if(int(row[1][-1]) >= 0 and int(row[1][-1]) <= 100 and int(row[0][-1]) >= 0 and int(row[0][-1]) <= 100):
    #        c_r = int(row[1][-1])
    #        c_g = int(row[0][-1])
    #        c_b = int(0)
    #        plt.plot(row[1],row[0], label='E', linewidth=2, color=(float(c_r/100), float(c_g/100), 0.0))
            # rE[0] is global and goes on the y axis


    for row in rE_prime:

        if((int(row[1][-1]),int(row[0][-1])) not in stable_locations):
            stable_locations.append((int(row[1][-1]),int(row[0][-1])))

        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            # float color must be between 0 and 1
            # trajectories outside 0 R
            plt.plot(row[1][0], row[0][0], marker='.', markersize = "10", color=(float(0), float(0), float(1))) # Plots Start but not the Ends
            plt.plot(row[1],row[0], label='E', linewidth=1, color=(float(0), float(0), float(1)))
        else:
            plt.plot(row[1],row[0], label='E', linewidth=1, color=(float(c_r/100), float(c_g/100), float(0.5)))

            plt.plot(row[1][0], row[0][0], marker='.', markersize = "10" , color=(float(c_r/100), float(c_g/100), float(0.5)))
            plt.plot(row[1][-1], row[0][-1], marker='*', markersize = "10" , color=(float(c_r/100), float(c_g/100), float(0.5)))

    plt.savefig("tra_reg_rgb" + str(RUN_ID) + "-" + str(random.randint(100, 999)) +".png")
    #plt.show()

#stable_points_space()



################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################



def stable_points_space_3d():

    fig = plt.figure(figsize=(30,30), dpi=200)
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Total Abundance")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Abundance', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])


# Abundance is captured like so:
    # [[1,2,3],[1,2,3],[1,2,3]... K]
    # and packed for each sample run into rAxS_prime
    # so should look like this
    #[
    #    [[1,2,3],[1,2,3],[1,2,3]... K],
    #    [[1,2,3],[1,2,3],[1,2,3]... K],
    #    [[1,2,3],[1,2,3],[1,2,3]... K],
    #    .
    #    .
    #    .
    #    SAMPLE_STEP
    #]

    index_A = 0

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_abundance = []
        row_abundance.append(0)
        decompose = rAx_prime[index_A]
        #print(decompose)
        run_length = len(decompose[0])
        #print(run_length)


        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_abundance.append(current_sum)

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_abundance[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #ax.scatter(row[1][0], row[0][0], row_abundance[0],s=10, marker='*', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], 0, marker='.', s=10, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_abundance[-1], marker='*', s=15, color=(float(c_r/100), float(c_g/100), float(0.5)))

    index_A +=1

    plt.savefig("3d_abundance_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    #plt.show()

#stable_points_space_3d()


def stable_points_space_3d_rotate():

    fig = plt.figure(figsize=(50,50))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Total Abundance")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Abundance', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])


    # Abundance is captured like so:
    # [[1,2,3],[1,2,3],[1,2,3]... K]
    # and packed for each sample run into rAxS_prime
    # so should look like this
    #[
    #    [[1,2,3],[1,2,3],[1,2,3]... K],
    #    [[1,2,3],[1,2,3],[1,2,3]... K],
    #    [[1,2,3],[1,2,3],[1,2,3]... K],
    #    .
    #    .
    #    .
    #    SAMPLE_STEP
    #]

    index_A = 0

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_abundance = []
        row_abundance.append(0)
        decompose = rAx_prime[index_A]
        #print(decompose)
        run_length = len(decompose[0])
        #print(run_length)


        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_abundance.append(current_sum)

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_abundance[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], 0, marker='.', s=10, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_abundance[-1], marker='*', s=15, color=(float(c_r/100), float(c_g/100), float(0.5)))

    index_A +=1

    #plt.savefig("3d_abundance_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

#stable_points_space_3d_rotate()

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
    s['rAx_prime']      = rAx_prime
    s['rAxR_prime']     = rAxR_prime
    s['time_prime']     = time_prime
    s['rE_prime']       = rE_prime
    s['rF_prime']       = rF_prime

    s['K']              = K
    s['R']              = R
    s['E']              = E
    s['start']          = start
    s['end']            = end
    s['step']           = step
    s['N']              = N
    s['OEn']            = OEn

    #s['a_t']                = a_t
    #s['a_l']                = a_l
    #s['a_g']                = a_g
    s['simulation_run']     = simulation_run
    s['RUN_ID']             = RUN_ID
    s['local_population_']  = local_population_

finally:
    s.close()

