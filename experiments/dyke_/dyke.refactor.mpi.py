import sys
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

if int(len(sys.argv)) != int(2):
    print("Args: shelve file name which contains all of > (K, R, P, E, start, end, step, EN, OE, LP_Z, RUN_ID)")
    print("e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN=2, OE=5, LP_Z = (10 - 100), RUN_ID : epoch")
    print("exit")
    sys.exit()

s = shelve.open(str(sys.argv[1]))

try:
    SAMPLE_SIZE = s['SAMPLE_SIZE']
    SAMPLE_STEP = s['SAMPLE_STEP']
    RUN_ID = s['RUN_ID']

    biotic_components_K = s['biotic_components_K']
    essential_range_R = s['essential_range_R']
    external_perturbation_rate_P = s['external_perturbation_rate_P']
    time_start = s['time_start']
    time_end = s['time_end']
    time_step = s['time_step']
    environment_components_N = s['environment_components_N']
    truncated_gaussian_ROUND = s['truncated_gaussian_ROUND']
    niche_width = s['niche_width']
    local_population_size = s['local_population_size']
    affects_w = s['affects_w']
    optimum_condition_u = s['optimum_condition_u']
    biotic_force_F = s['biotic_force_F']

    exp_name = s['exp_name']
    data_directory = s['data_directory']
    shelve_file = s['shelve_file']

finally:
    s.close()

S_STEP = SAMPLE_STEP
K = biotic_components_K
R = essential_range_R
P = external_perturbation_rate_P
start = time_start
end = time_end
step = time_step
N = environment_components_N
E = [0, 0]
F = biotic_force_F

ROUND = truncated_gaussian_ROUND

OEn = niche_width
OE = [OEn for _ in range(K)]

local_population_ = local_population_size
local_population_index = []
local_population_size = int(local_population_/100 * K)

uniq_k = []
for x in range(local_population_size):
    one = random.randint(0,K-1)
    while one in uniq_k:
        one = random.randint(0,K-1)
    uniq_k.append(one)
    local_population_index.append(one)


local_population_index.sort()

w = affects_w
u = optimum_condition_u

w = [[0.6037119304464391, -0.9889789963401783, -0.29484339996405895, 0.364933779603988, 0.0775549915728031, -0.5580454891206899, 0.11289892783260513, 0.15560096870687978, -0.29897915018120313, -0.01596015326011191, 0.31616464757050267, -0.8949950896708554, -0.8635823255212485, 0.3425326625828171, -0.4742372122247336, -0.48867918454090864, 0.48043588476553767, -0.760410057029262, -0.6646405760277712, 0.8407561462658837, -0.28194705759505423, -0.47465080030926576, -0.4099085785830754, 0.9890104400713089, -0.7753485721536115, 0.059295909818493975, -0.44066917983204323, 0.26468440374931035, 0.865083764178955, -0.5768419729493097, -0.5753904807105319, -0.875183980364524, -0.2886232189818301, 0.1453132194576119, -0.07854162342300142, -0.7931700432501925, -0.2636277485651546, -0.513063287150737, 0.5431629548838461, -0.1414061761578307, 0.09091890867306662, 0.5731239838804805, 0.4015113299254589, -0.897834425630039, 0.6815699293088158, -0.04412309184192642, 0.5337184094195873, 0.0009570117759403196, 0.9527757520727242, 0.9847963860569273, 0.8535623930616967, 0.6602490223279283, 0.6094154432249821, 0.9912986651145881, -0.7058415095415864, 0.04163235694546974, 0.3572733519195479, -0.2938892591978368, -0.9455533780806107, 0.2151273038344339, 0.7896590907156757, -0.16206913155099034, 0.8896460138395927, -0.1608253970843856, -0.1943344406890546, -0.4504035056920117, 0.06895260964524041, 0.9661839300258792, 0.5090686511942188, 0.49072669171365013, 0.20543671404088704, 0.12322015615195459, 0.6668940558915815, -0.8625724730201931, 0.7851717005928216, -0.39451535370609814, -0.23401296802072125, 0.1591918253929967, -0.14163768013658973, -0.0876951470864189, -0.9684866275673447, -0.7104854461481509, -0.46370128942652333, -0.04261219481265521, 0.34418141179925543, 0.35598764870718647, -0.9486849443183756, 0.32447567584022274, 0.3467546016304024, -0.7160578344581641, -0.5257070301706246, 0.1932142677099704, 0.37526078981738475, -0.8570082699144232, 0.22045072452413628, 0.8036137207716259, -0.7991339139426676, 0.8237234128242117, -0.7936307473019131, 0.3415459423677365], [-0.6772300105464732, -0.8760694091426917, -0.5963903101258312, 0.7791544065818952, 0.05364230028168415, 0.8177403988944534, 0.8876103692454418, -0.6014953672883911, 0.0040660721845671155, -0.7700312921709054, 0.2749823043585349, -0.4790660105407738, 0.25861507950321005, 0.9079601286077468, 0.14692686083113093, 0.03419378245425664, 0.10396076110015118, 0.038980896667212495, 0.13071137250641152, 0.98537508499391, 0.4162582367024008, 0.9444649713195066, -0.8255291184117568, 0.4847639967669648, 0.5222062852991949, 0.14170089787321682, -0.46807069303830406, -0.05283852529299726, 0.4664699873977469, 0.9882325638497489, 0.8873686905866611, 0.9656604259513457, -0.47263307121569387, 0.8430954422549837, -0.8296528302817507, 0.8372885365380456, 0.13865367029516684, -0.44723865521808626, -0.9220753216783153, -0.3588767897549612, -0.10964463671617275, -0.1987452395269127, 0.856127216733424, -0.8305021615901425, -0.08866289548990158, -0.014006988595108849, 0.7977534488802827, -0.1692953475990513, -0.0015953029547850495, -0.20201494420343735, -0.9928238873128497, 0.35741971342192946, -0.20233783134514738, -0.09345313656175636, 0.3038524498259618, 0.48092446767517294, -0.227852956259754, -0.4323264303060905, -0.26191045107166167, 0.826660957685659, 0.4894683511964899, -0.8703329123427694, -0.8781357944843176, 0.04607679150858601, -0.291270050860438, -0.788399069753627, -0.48679101649720224, 0.6924066733440919, -0.5480758953093281, 0.1907023585086911, 0.6911389861982482, 0.40186494457960453, 0.6380250530984299, 0.40049804531211275, 0.8139163542369503, 0.783333939882326, -0.04590939481533418, 0.33937968494523707, -0.43569592268018376, -0.41439350794144336, -0.6134954756057034, 0.3406775832829827, -0.15489107889501552, 0.3929436314404182, -0.15209714289370124, 0.2532014215236027, -0.3546215903122367, 0.9262752425169627, 0.059225160529298604, -0.06952470223165652, 0.8539019334906339, -0.12274957744400239, -0.7135515740341996, 0.3504215976195968, 0.27452853668650623, -0.5017098896802097, 0.35084841173881354, -0.5715938790597674, 0.8951212681283327, -0.2958294137343982]]
u = [[23.721279453271016, 47.047876088197526, 53.51752340497219, 75.22875314386675, 38.65686647546275, 26.79825817954613, 78.77880861842976, 68.88683520043398, 71.61663725544373, 64.5332702128044, 22.478320111985706, 96.72046585552133, 44.00805947168222, 77.34899086597031, 59.1011949116247, 47.88262420014876, 60.689965502149015, 65.31195674315407, 53.35037343816966, 14.306712257351052, 84.20703827057751, 68.70715774932455, 60.28955886937661, 11.073636384783148, 26.14559225934193, 64.09421016349664, 39.57727886783751, 25.590543450334536, 0.4121493398364984, 92.88119928952348, 31.484460367121102, 0.10685735701121413, 26.589445045692973, 16.674600181226452, 93.57239628225585, 60.833164657416184, 77.86396606214699, 19.24738816906453, 66.60547231230326, 73.58005060998619, 75.42681024498037, 8.599519877124484, 55.828379075651334, 84.17896731914888, 32.15675729637757, 94.28926184672576, 80.32204470067906, 71.97086764118824, 82.77764322905189, 16.095287732421415, 42.09553060517109, 73.75276402725062, 93.9142591129815, 77.14845583685052, 64.02529952059646, 93.03317634750717, 9.350643682384241, 51.50809916288732, 95.90149776699108, 57.233501265462635, 11.128126044250985, 43.08007916002868, 71.32557020389089, 62.06461923009311, 57.15271195310543, 6.035417764721551, 54.792846370948425, 96.47140825696256, 15.578520806351193, 83.67485939920203, 67.89223790752665, 70.63648340502037, 17.59998959597843, 12.234260764167892, 67.03729671349156, 20.845006587113414, 37.42770879432909, 43.67298649733636, 44.53797660611335, 13.591866348979742, 46.816788052074955, 59.94604061035154, 25.253715222233886, 93.33471407972824, 60.78839024536254, 77.73794549355722, 90.36663035884807, 61.92853257847178, 33.23124800222396, 97.18921145873252, 60.60943833852326, 80.31272271809976, 65.77066031481745, 67.88585127924225, 1.9204548373100039, 53.50801173203745, 18.47861385347702, 65.76890193857695, 81.27169748755703, 85.1543471082693], [49.72596191463289, 82.47171943590327, 40.57892655117874, 79.82044675569297, 82.38792289965812, 37.851846852508245, 43.41288534133143, 27.070412613010674, 52.69631231091836, 3.381831592670992, 21.499624142670505, 73.06227148593346, 31.539500071585547, 72.7607297623941, 90.75872811545295, 34.750783830527844, 0.359239314110682, 47.555397598980434, 34.62304450951799, 26.379637190761883, 58.6385279181903, 18.21744635064896, 92.23744211120474, 57.59623858143167, 45.56668800884479, 55.19345539753123, 49.79758384502697, 64.55104512301237, 72.93523558565167, 95.46874045036783, 30.0208420153651, 46.930160779945496, 13.410917769177576, 98.23725790283221, 82.25485073795205, 45.88716207467396, 80.72192023087229, 50.95617824019105, 2.1440681873308387, 28.468845989440794, 76.44014716084388, 94.17661008651822, 67.38042071534947, 15.854433827976955, 61.102867952277926, 78.5485596981821, 70.85876786366109, 53.13503856097746, 98.95651702260567, 83.01118534675334, 65.32638931646477, 89.15893243223226, 53.6927400554229, 63.20597675136851, 8.796668198413137, 82.85610880635619, 76.24257935686107, 57.49716757183675, 5.4571367131231, 24.337036025455816, 27.363318995169838, 76.69481067844457, 12.079727345281565, 2.1920168457766565, 93.15613342073192, 47.99081951796956, 77.63435619729125, 32.540590966953786, 55.84659769119662, 0.2634949266941078, 49.191901917583436, 78.14274947977033, 90.21790346735423, 46.19693373995533, 98.64646136399058, 79.61627443746342, 34.72193246428582, 18.732709509437573, 36.131249540536146, 26.534599297233697, 81.5244124870408, 61.44163455765476, 94.21407279218121, 51.86159537037434, 41.51587910936653, 20.193804607196665, 71.50665280709478, 30.91283943546682, 84.19691593025156, 94.28964994005796, 11.32682153580159, 10.449828543371009, 66.83727630057034, 97.8031295261222, 80.52856095374936, 51.253788797679746, 28.454303939634183, 9.46271787709242, 94.74939900058112, 59.09552824656377]]


print(w)
print(u)

#Abundance values over time
rAx = [[] for x in range(K)]
#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]

rNumberAlive = [[] for x in range(K)]

#for _ in range(K):
#    rNumberAlive[_].append(0)

alpha = [[] for _ in range(K)]

# Fix applied - correctly generating alphas for global and local
for _ in range(K):
    al = []
    for ai in range(N):
        al.append(round((math.e) ** ((-1) * (((abs((E[ai])-u[ai][_])) ** 2) / (2*(OE[_]**2)))),ROUND))

        new_alpha = 0
        if _ in local_population_index:
            new_alpha = np.prod(al) # If local population (product of abundances) (both local and global affect the local one)
        else:
            new_alpha = al[0]       # Else take the first one as Eg

        alpha[_].append(new_alpha)

        rAx[_].append(new_alpha)
        rAxR[_].append(new_alpha * R)

        if(new_alpha > 0):
            rNumberAlive[_].append(1)
        else:
            rNumberAlive[_].append(0)

print("alpha: ",alpha)
print("Number Alive : ", rNumberAlive)

rF = [[] for _ in range(N)]         #Biotic Force Values
rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

initfSUM = [0 for _ in range(N)]
for _ in range(K):
    initfSUM[0] = initfSUM[0] + (alpha[_][-1] * w[0][_])

rF[0].append(initfSUM[0])

for _ in range(K):
    if( _ in local_population_index):
        initfSUM[1] = initfSUM[1] + (alpha[_][-1] * w[1][_])

rF[1].append(initfSUM[1])


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

        new_alpha = 0
        if( _ in local_population_index):
            new_alpha = np.prod(al)
        else:
            new_alpha = al[0]

        k1 = new_alpha
        k2 = k1 + (k1 * step/2)
        k3 = k1 + (k2 * step/2)
        k4 = k1 + (k3 * step)
        yt = alpha[_][-1] + (((k1 + (2*k2) + (2*k3) + k4)/6) - alpha[_][-1]) * step

        alpha[_].append(yt)

        rAx[_].append(yt)
        rAxR[_].append(yt * R)

        if(yt > 0):
            rNumberAlive[_].append(1)
        else:
            rNumberAlive[_].append(0)

    for _ in range(K):
        fSUM[0] = fSUM[0] + (alpha[_][-1] * w[0][_])

    F[0] = fSUM[0]
    newE = E[0] + (((0 + F[0]) * 1))
    E[0] = E[0] + ((newE-E[0]) * temperature_time_scale)
    rF[0].append(F[0])
    rE[0].append(E[0])

    # ============ END EG ==============

    for _ in range(K):
        if( _ in local_population_index):
            fSUM[1] = fSUM[1] + (alpha[_][-1] * w[1][_])

    F[1] = fSUM[1]
    newE = E[1] + (((0 + F[1]) * 1))
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
    rNumberAlive_prime  =[]


    # First Set of calculations have occurred during initilization so appending time 0


    # sampling
    for Eg_temp in np.arange(1,100,S_STEP):
        for El_temp in np.arange(1,100,S_STEP):
            print("Init : ", Eg_temp, El_temp)
            simulation_run.append((Eg_temp,El_temp))
            time.append(0)
            # xtime should should start from one timestep + 0
            post_init_start = start + step
            for xtime in np.arange (post_init_start, end, step):
                update(step)
                time.append(xtime)

            E_prime.append([Eg_temp, El_temp])
            F_prime.append(F)
            for _ in range(K):
                del alpha[_][-1]
            alpha_prime.append(alpha)
            rF_prime.append(rF)
            rE_prime.append(rE)
            for _ in range(K):
                del rAx[_][-1]
            rAx_prime.append(rAx)
            for _ in range(K):
                del rAxR[_][-1]
            rAxR_prime.append(rAxR)
            for _ in range(K):
                del rNumberAlive[_][-1]
            rNumberAlive_prime.append(rNumberAlive)
            time_prime.append(time)
            biotic_force_prime.append(biotic_force)
            temperatures_prime.append(temperatures)

            #fig = plt.figure(figsize=(50,50))
            #ax = fig.add_subplot(111, projection='3d', adjustable='box')
            #ax.set_title(label = "EL/EG with Number Alive Single")
            #ax.set_xlabel('X - EL', fontsize=10)
            #ax.set_ylabel('Y - EG', fontsize=10)
            #ax.set_zlabel('Z - Total Alive', fontsize=10)

            #c_r = int(rE[1][-1])
            #c_g = int(rE[0][-1])

            #print(len(rE[0]))
            #print(len(rE[1]))
            #print(len(rNumberAlive))
            #print(len(rNumberAlive[0]))
            #sum_alives = []
            #for x in range (len(rNumberAlive[0])):
            #    current_sum = 0
            #    for _ in range(K):
            #        current_sum += rNumberAlive[_][x]
            #    sum_alives.append(current_sum)

            #print(len(sum_alives))
            #if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            #    ax.scatter(rE[1][0], rE[0][0], sum_alives[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #    plt.plot(rE[1],rE[0], sum_alives , color=(float(0), float(0), float(1)))
            #else:
            #    plt.plot(rE[1],rE[0],sum_alives, color=(float(c_r/100), float(c_g/100), float(0.5)))
            #    ax.scatter(rE[1][0], rE[0][0], 0, marker='.', s=10, color=(float(c_r/100), float(c_g/100), float(0.5)))
            #    ax.scatter(rE[1][-1], rE[0][-1], sum_alives[-1], marker='*', s=15, color=(float(c_r/100), float(c_g/100), float(0.5)))

            #plt.show()

            ###########################################################################################################################
            ###########################################################################################################################
            ###########################################################################################################################
            ######################################### RE INIT #########################################################################
            E = [Eg_temp,El_temp]
            F       = [0 for _ in range(N)]
            alpha = [[] for _ in range(K)]
            rAx = [[] for x in range(K)]
            rAxR = [[] for x in range(K)]
            rNumberAlive = [[] for x in range(K)]

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
                    rAx[_].append(new_alpha)
                    rAxR[_].append(new_alpha * R)
                    if(new_alpha > 0):
                        rNumberAlive[_].append(1)
                    else:
                        rNumberAlive[_].append(0)

            rF = [[] for _ in range(N)]         #Biotic Force Values
            rE = [[] for _ in range(N)]         #A blank list for each Environment Variable

            initfSUM = [0 for _ in range(N)]
            for _ in range(K):
                initfSUM[0] = initfSUM[0] + (alpha[_][-1] * w[0][_])

            rF[0].append(initfSUM[0])

            for _ in range(K):
                if( _ in local_population_index):
                    initfSUM[1] = initfSUM[1] + (alpha[_][-1] * w[1][_])

            rF[1].append(initfSUM[1])

            for _ in range(N):
                rE[_].append(E[_])                 #Input the Start Temperatures




            time = []
            biotic_force = [[] for _ in range(K)]
            temperatures = []

            ###########################################################################################################################
            ###########################################################################################################################
            ###########################################################################################################################
            ######################################### END RE INIT #####################################################################


def plot_alphas():

    for x in np.arange (-50, R+50, step):
        temperatures.append(x)

    plt.figure(figsize=(30,30))
    plt.title('Biotic Force over Temperature', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    super_biotic_force = []

    for each_env_var in range(N):
        biotic_force = [[] for _ in range(K)]
        for y in range(K):
            for x in np.arange (-50, R+50, step):
                biotic_force[y].append((math.e) ** ((-1) * (((abs(x-u[each_env_var][y])) ** 2) / (2*(OE[y]**2)))) * w[each_env_var][y])


        for _ in range(K):
            plt.plot(temperatures,biotic_force[_])


        plt.plot(temperatures,np.sum((np.array(biotic_force, dtype=float)), axis=0), lw=4)
        super_biotic_force.append(np.sum((np.array(biotic_force, dtype=float)), axis=0))

    sum = []
    #for _ in range(time_end/time_step):
    #    sum.append(super_biotic_force[0][_] + super_biotic_force[1][_])
    #plt.plot(temperatures, sum, lw=10)
    #plt.show()
#plot_alphas()


def stable_points_space():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Regions', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)


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

    #plt.savefig("tra_reg_rgb" + str(RUN_ID) + "-" + str(random.randint(100, 999)) +".png")
    plt.show()

stable_points_space()


def stable_points_space_heat():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Heat Map of Global and Local Species Distribution', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    # Species surviving
    index_global=0

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for optimum_condition in u[0]: # Global
        plot_x=[]
        plot_y=[]
        if(index_global not in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1

    index_global=0
    for optimum_condition in u[1]: # Local
        plot_x=[]
        plot_y=[]
        if(index_global in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1

    plt.imshow(heatmap, cmap='hot', interpolation='nearest')

    plt.show()

stable_points_space_heat()

def stable_points_space_heat_global():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Heat Map of Global Species Distribution', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    # Species surviving
    index_global=0

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for optimum_condition in u[0]: # Global
        plot_x=[]
        plot_y=[]
        if(index_global not in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1

    plt.imshow(heatmap, cmap='hot', interpolation='nearest')

    plt.show()

#stable_points_space_heat_global()

def stable_points_space_heat_local():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Heat Map of Local Species Distribution', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)
    # Species surviving
    index_global=0

    heatmap = [[0 for _ in range(R+7)] for _ in range(R+7)]

    for optimum_condition in u[1]: # Local
        plot_x=[]
        plot_y=[]
        if(index_global in local_population_index):
            for x in np.arange(0,R,1):
                for y in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            plot_x=[]
            plot_y=[]
            for y in np.arange(0,R,1):
                for x in np.arange(optimum_condition-OEn,optimum_condition+OEn,1):
                    plot_x.append(x)
                    plot_y.append(y)
            local_index = 0
            for location in plot_x:
                heatmap[int(plot_x[local_index])][int(plot_y[local_index])] += 1
                local_index+=1
            #plt.scatter(plot_x, plot_y, marker='o', alpha=0.03, color="grey", edgecolors='none')

        index_global +=1

    plt.imshow(heatmap, cmap='hot', interpolation='nearest')

    plt.show()

#stable_points_space_heat_local()

def stable_points_space_final_abundance_rgb():

    stable_locations = []
    plt.figure(figsize=(30,30), dpi=200)
    plt.title('Final Abundance Map', fontsize=40)
    plt.xlabel('EL', fontsize=40)
    plt.ylabel('EG', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-20, R+20)
    plt.xlim(-20, R+20)

    index_A = 0
    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_abundance = []
        decompose = rAx_prime[index_A]
        run_length = len(decompose[0])

        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_abundance.append(current_sum)

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            plt.scatter(row[1][0], row[0][0],s=10, marker='.', color=(float(row_abundance[0]), float(1), float(1)))
            #plt.plot(row[1],row[0], row_abundance , color=(float(0), float(0), float(1)))
        else:
            #plt.plot(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)))
            #plt.scatter(row[1][0], row[0][0], row_abundance[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            plt.scatter(row[1][-1], row[0][-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(row_abundance[-1])))

        index_A +=1
    plt.show()

stable_points_space_final_abundance_rgb()

mpl.use('macosx') #for the 3D magic

def stable_points_space_3d_rotate():

    fig = plt.figure(figsize=(300,300))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Total Abundance")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Abundance', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])

    index_A = 0

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        #print(len(row[1]))
        #print(len(row[0]))
        #print(len(rAx_prime[index_A]))


        row_abundance = []
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
            ax.scatter(row[1][0], row[0][0], row_abundance[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_abundance[-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))

        index_A +=1

    #plt.savefig("3d_abundance_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

stable_points_space_3d_rotate()


def stable_points_space_3d_rotate_number_alive():

    fig = plt.figure(figsize=(300,300))
    ax = fig.add_subplot(111, projection='3d', adjustable='box')
    ax.set_title(label = "EL/EG with Number Alive")
    ax.set_xlabel('X - EL', fontsize=10)
    ax.set_ylabel('Y - EG', fontsize=10)
    ax.set_zlabel('Z - Total Alive', fontsize=10)
    ax.set_xlim([-10,100])
    ax.set_ylim([-10,100])
    #ax.set_zlim([0,50])

    index_A = 0

    #print(rNumberAlive_prime)

    for row in rE_prime:
        c_r = int(row[1][-1])
        c_g = int(row[0][-1])

        row_alive = []

        decompose = rNumberAlive_prime[index_A]
        run_length = len(decompose[0])

        for x in range(run_length):
            current_sum = 0
            for item in decompose:
                current_sum += item[x]
            row_alive.append(current_sum)

        #print(len(row_abundance))
        #print(len(row[1]))
        #print(len(row[1]))

        if(c_r < 0 or c_g < 0 or c_r > 100 or c_g > 100 ):
            ax.scatter(row[1][0], row[0][0], row_alive[0],s=10, marker='.', color=(float(0), float(0), float(1)))
            #plt.plot(row[1][0], row[0][0], row_abundance[0], marker='x')
            plt.plot(row[1],row[0], row_alive , color=(float(0), float(0), float(1)))
        else:
            #ax.scatter(row[1],row[0],row_abundance, color=(float(c_r/100), float(c_g/100), float(0.5)), s=1)
            plt.plot(row[1],row[0],row_alive, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][0], row[0][0], row_alive[0], marker='.', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))
            ax.scatter(row[1][-1], row[0][-1], row_alive[-1], marker='*', s=20, color=(float(c_r/100), float(c_g/100), float(0.5)))

        index_A +=1

    #plt.savefig("3d_alive_" + str(RUN_ID)  + "-" +  str(random.randint(100, 999)) + ".png")
    plt.show()

stable_points_space_3d_rotate_number_alive()