import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
from matplotlib.gridspec import GridSpec

import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
#from numba import jit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
plt.rcParams["font.family"] = "Times New Roman"


# Generating ALL Parameters
SAMPLE_SIZE = 1
SAMPLE_STEP = 1
RUN_ID = int(time.time())

SPECIES_K   = 100                  # ----------- Number of Biotic Components
RANGE_R     = 100                  # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 1                   # ----------- Time Step3
ENV_VARS    = 1                     # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0
ENV_START=[10]
omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

system_state = np.zeros(SPECIES_K+ENV_VARS)

# system_state[0-99]  = species
# system_state[100]   = environment
# system_state[101]   = biotic_force
# system_state[102]   = perturbing force
# system_state[103]   = exponential temp runoff
global biotic
biotic = []
global perturb
perturb = []
global exp_temp
exp_temp = []
global add
add = 0

Eg = ENV_START[0]

for s_i in range(SPECIES_K):
    a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

    if a_star < ALIVE_THRESHOLD:
        a_star = 0

    system_state[s_i] = a_star

for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START[_]

def rates_of_change_system_state(system_state):

    rate_of_change = system_state.copy()

    Eg = system_state[SPECIES_K+0]
    #ABUNDANCE
    for s_i in range(SPECIES_K):

        a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))

        rate_of_change[s_i] =  a_star - system_state[s_i]

    biotic_force_FG = 0
    #BIOTIC FORCE
    for s_i in range(SPECIES_K):
        biotic_force_FG += (system_state[s_i] * omega[0][s_i])
    rate_of_change[SPECIES_K+0] = (biotic_force_FG)

    if(add == 1):
        #BIOTIC
        biotic.append(biotic_force_FG * 10)
        rate_of_change[SPECIES_K+0] = (biotic_force_FG + perturb[-1])
        #TEMP
        exp_temp.append(exp_temp[-1] + (exp_temp[-1] + perturb[-1])/20)
        #PERTURBING
        perturb.append(perturb[-1]+0.2)

    return(rate_of_change)

if __name__ == '__main__':


    # system_state[0-99]  = species
    # system_state[100]   = environment
    # system_state[101]   = biotic_force
    # system_state[102]   = perturbing force
    # system_state[103]   = exponential temp runoff

    for x in range(1):
        omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        # ======================================================================
        # 39
        omega = [[-0.40382416079337524, 0.06387871687541957, 0.7171761889407076, -0.8824867928095121, 0.609701346229401, -0.9778612528385153, -0.6904390349071134, -0.7293413088469507, 0.8718318987555833, -0.11805319404462411, -0.3222583844152793, 0.26008251337899035, 0.6652738211270313, -0.7439556291997464, -0.561165123492809, 0.8107326144697891, -0.29117921839816896, -0.09250315285609223, -0.04295473401061467, -0.4210937587465704, -0.2177254843097034, -0.044209883925610205, -0.4515411507451814, 0.01836111172374233, 0.21854662961108384, 0.024021265526827484, -0.9021880519407839, 0.8336320137535931, -0.12932861424360764, -0.2406540490601119, 0.9056198299036315, -0.3111922610493725, 0.9932853941430557, -0.0028959529243335336, -0.9783497699858072, -0.645599573001375, 0.326735869753366, 0.2473307843837822, -0.8349734263711832, 0.8753052137857333, -0.1429988993103959, -0.5889554435241113, 0.7787334370668857, 0.2755842899474743, 0.7240351333602393, 0.23217717644810287, -0.6644413576212658, -0.19136483822052108, -0.5143458135261845, 0.3895213411272149, -0.1645177935599258, -0.3993998474371723, -0.7377666841605925, 0.6733658024513127, 0.4612162049117312, -0.8146815813280639, -0.670965586335025, -0.5093531941189002, 0.6191083735282072, -0.2739295032486042, -0.5211315755757249, -0.8299901864561445, -0.5795596343420406, -0.3293599032960408, 0.1769052749150155, 0.3883385408675235, 0.6682779918709842, 0.7715538954023624, -0.4107023013267388, -0.06576985249889189, -0.13009081522919086, 0.9695666038466064, -0.12421010213951167, -0.12206838826739119, 0.04159413271339707, 0.01366272161394022, -0.64156947448589, -0.35497037061038905, 0.5791602152488247, -0.9843331717596961, -0.9228105122483896, 0.1348397221386508, -0.692064757893107, 0.6091106291710222, 0.008692211030396635, -0.6236024692766666, 0.9984311866444902, -0.551444385196312, -0.07308760911476164, -0.6965614624760834, 0.12800805134247129, 0.8129015896054292, 0.737294771197603, -0.5753767306709692, 0.5803218786218514, 0.4422015191031661, -0.0997702044715838, -0.7505586092626417, 0.817111842989015, -0.23066051842257762]]
        mu = [[57.70150907906662, 64.95894538234288, 53.88534521282694, 92.54478018625232, 77.58414205145534, 21.67956140667152, 20.663109763418653, 1.0552597756285809, 91.38038123664305, 75.31082151307149, 72.0217722019113, 43.75199886287102, 11.233724196004191, 89.16756105624651, 90.19828138052651, 52.86990946240447, 54.2413254445436, 75.13697352248502, 98.12068008206451, 64.55163539014791, 66.13308435447837, 72.59910438585881, 47.075273400181125, 40.828401826098684, 40.52387810523522, 2.533052986312112, 16.67812283945701, 20.319526222324036, 45.21279623969545, 93.15634683271655, 13.488847865121734, 24.434080323122387, 81.58624243610933, 46.935186765276505, 58.019708238631985, 85.99494640688924, 65.79213485603886, 79.50489548207199, 56.09890292752252, 79.56587422354391, 93.10427753207846, 50.49747812732414, 31.91489111321434, 24.07638414095652, 50.52171624210803, 54.86912249758515, 66.5277652476428, 10.413103661751988, 64.26899310814235, 8.471240846685024, 49.691852799676695, 71.8505894583984, 65.73236135367522, 96.94189783466963, 50.32638050564516, 20.975647583224124, 62.98387776922484, 27.643824332970823, 95.18096250614423, 74.99739028713992, 58.699666385612346, 60.18474578277214, 99.60902173491964, 94.10258650011274, 14.45760677324186, 16.56316599037193, 79.1756194185021, 69.31905962638005, 74.78938514478897, 70.43603863435395, 43.77075308908108, 16.59872955445959, 78.6913245570703, 93.63391331741332, 25.375275245211615, 12.242956789118065, 62.91091600521619, 2.1647808864662776, 17.59895352715609, 17.796251411202967, 20.424276493307925, 47.92090295614584, 23.967810729587157, 71.81090549506624, 35.626995477250304, 50.75762668559681, 24.364775064805855, 7.747601519641201, 4.603022727301886, 92.7582813329584, 79.64283032203878, 39.07341098148338, 76.25017342615332, 35.24521767497608, 30.369229065577375, 82.27974647973346, 44.2836989830251, 48.25819736008733, 30.690072664351153, 3.7007136722621903]]
        # ======================================================================
        # 64
        #omega = [[-0.9927479585284591, 0.4051662451289728, -0.3991513631186596, 0.988227260134791, 0.02717414453370015, -0.7978452099867359, -0.6517727574737826, 0.1031684585525301, 0.8303285224358066, 0.7982319388441845, 0.4176436306889413, 0.28919010826746194, 0.17379339913813885, -0.8006995782041539, 0.7988324215105351, 0.5639747767404775, -0.18337101655189603, 0.41573385248738925, -0.3265318912186168, -0.8019454168036535, 0.27116883483087295, 0.9312280517007812, -0.8672210419552067, -0.9061091180169782, -0.30891828327247906, 0.3485172406111503, 0.8668344586195265, 0.6364317190848847, -0.5982091566484828, -0.0612532502591534, 0.8597241964428721, -0.6522145988673882, 0.20411409884756138, -0.7478520586370627, 0.9052655158333613, -0.5492509876203651, -0.8685223147139627, -0.32460542741277676, 0.26052515149039057, 0.005302343889805661, 0.7203004804793172, -0.09199414893928792, -0.5765853359027504, 0.02749704998792457, 0.049110107417962734, -0.37305815110897855, -0.6617943166020137, 0.02921049520073149, -0.09310071216197935, 0.18015633863523317, 0.15808813008590028, 0.4738221452810467, 0.8976843632519982, -0.19598626106494121, -0.4259566495569158, 0.03280109021895705, -0.09475667773819096, -0.30331200449309836, -0.7565128774740024, -0.07524825944923585, -0.03058152312416218, -0.7293251146353206, 0.48163006871219616, -0.5467965910976684, -0.13565519833239437, -0.21950059214215245, -0.8558154779505052, -0.590269873242643, -0.03845533402680679, 0.34093980016018777, 0.46340540146788367, -0.5740989729841341, 0.9086765251523614, 0.27121395514098423, 0.35987346231938355, -0.17382170780371764, -0.971122926456232, -0.12641225777721776, -0.0052092933248031326, 0.1842996340127876, -0.8436748964012459, 0.24743034614908943, -0.7996607925590209, -0.004665257127270284, -0.17418890263514908, -0.5488680661263945, 0.4757360491789999, -0.7366587821115762, 0.2829999474920497, -0.5301925267039038, -0.697313497820713, 0.2712725101960294, 0.5291964531778737, 0.9781210404687326, -0.08577286822267816, 0.5298901715027666, -0.26036629635641195, 0.03594746232816437, 0.41608575185500674, -0.5212984818158624]]
        #mu = [[68.03556389632452, 50.038953437338385, 17.405796053703227, 36.77431882772457, 43.53154508313134, 84.4011606159107, 31.359177382918634, 35.28775735156168, 17.82728927445446, 11.770777348307448, 19.53153237116322, 7.967051971640037, 77.06804927473489, 88.62649326504483, 72.15938367236888, 93.19668002096499, 78.9556727398785, 75.42576116348141, 6.040267219929973, 2.3728725115468197, 40.93957606117944, 31.09394674788486, 8.398152141104065, 57.57793730048263, 69.32166231955854, 3.602400955564, 17.500653956709257, 20.382165103459283, 80.49890083231828, 27.005669141075096, 99.95231232560388, 52.05416252355525, 64.03942700458643, 45.230362331813076, 45.074528044239436, 61.862338951987795, 98.51032396339996, 1.0599609826256295, 41.34845713591893, 3.9065022535100358, 40.69350968956685, 40.259230733266605, 83.80775296575214, 80.75739738578149, 2.631997942181874, 31.092316877336714, 55.71069782236234, 51.00669634528112, 8.883281597198822, 76.74339611912332, 29.7336367995404, 77.18414480789134, 95.75052556059708, 83.9513740773295, 30.16365226964186, 7.891212852454443, 50.22132027124747, 7.132815078581256, 93.21996607471425, 83.61411120222691, 58.06284365388499, 24.923781581076334, 61.35554350011824, 87.71060459677015, 2.468853870999277, 56.576496421045455, 28.945465866211407, 91.55857537894776, 55.229330990680246, 25.593862543615785, 47.77170876005362, 79.42080198711807, 29.843062702707513, 96.42765664498039, 32.365683710474016, 77.85626761183121, 49.123573655634466, 94.70021826496559, 14.441526840021545, 93.10500221788489, 20.201831163662078, 93.41215009015158, 93.4576007020496, 58.58347642102068, 90.93692022639992, 7.389895072061858, 31.3834450071083, 22.87984355432321, 20.379034232082216, 76.52955699967328, 25.478131953497083, 31.416558402319474, 31.70321279728049, 8.920234889440826, 4.136721897613649, 32.69675762807477, 15.997178808314395, 63.16352457348623, 50.219476856504585, 50.88301414065805]]
        # ======================================================================
        # 92
        #omega = [[0.6313967534802143, 0.8153593666944567, -0.6490534390931613, -0.5928271457380083, -0.6243231017738855, -0.3777583075799087, -0.35039455929416796, 0.9802572228875512, 0.40783495164632155, 0.218750884455708, 0.11069609037927108, 0.0669788312044779, -0.16962059475257107, -0.6893816627481357, 0.51218813837943, 0.7113111052799657, -0.012294751887598432, -0.4090898582750824, -0.14935258255008033, 0.03764607001770326, 0.08188955378508478, 0.5882649476867761, 0.35858367063763663, 0.8042120875635235, -0.21503324845294536, -0.5584452849725501, 0.560339872884593, 0.3087223067707958, 0.8571163524655334, -0.17072187729146648, 0.9333208641480941, 0.31707335604077014, 0.5477758781479165, -0.5907324125727575, -0.08607489788889389, 0.5964351300747981, -0.7092214684854583, 0.31860311277333664, -0.9749579989556998, -0.11543137759262345, 0.5826539436957174, -0.7409723723964414, -0.1754243000233855, -0.3000412727963173, 0.8430448057369806, 0.41806878328078945, 0.9955552667287104, 0.5617623037336321, -0.47382284740491953, -0.8400582187989112, 0.292825504434328, -0.8702874383153927, 0.38877111521414065, -0.48559187626505884, -0.18765117465662473, 0.4720281583851982, 0.7593086184890838, -0.43997132994659194, -0.6767707553444475, 0.597453539324099, -0.017529613122354126, -0.5114643856268168, -0.09586860341168357, 0.9785737019141347, -0.2162587543204535, 0.20704689694882417, 0.6991197760243142, 0.5479365241266152, -0.6576927315478303, -0.5716660588619173, -0.5886461267618175, 0.200718354066606, -0.6532850438767099, -0.6850373917028609, 0.703522531631988, 0.03135169508060698, -0.606150669766615, 0.5432997303728166, -0.5142905027813978, 0.48702931313526654, 0.734399065023208, -0.15258723265651253, -0.39461125161731636, -0.965215521024599, 0.5951408310203881, -0.24097252489237309, -0.019371190353112633, -0.8493888363934268, -0.17277940868519703, 0.5155970583616465, -0.36046301880926657, -0.019596435264123357, 0.6073102021608086, 0.12712690628556933, -0.9813499229228928, -0.6432964428506578, -0.733844989535769, 0.5249732042939805, 0.8962060131426575, -0.7766273067833269]]
        #mu = [[37.94906219174737, 41.51043119094561, 35.12300426349807, 50.21102451511516, 82.09333961269807, 9.946006478581037, 62.6318415326826, 76.18326306976635, 17.580260808027994, 36.83307385613281, 30.972290965212345, 15.855839150672468, 57.30267470463439, 94.71518376771657, 27.23466068486864, 12.291855674999553, 85.87007009580863, 30.15207395135532, 38.48405759977573, 10.2276218631829, 2.6049710071815113, 49.121699869410584, 92.47050141783079, 77.30702078738885, 6.779930144003254, 52.6120417084988, 75.35181106485214, 83.8503150516, 94.84632004357884, 46.80065177602744, 57.90624188009831, 71.29173774769014, 63.55610996775892, 8.106729586350337, 66.23770651906774, 11.31361777355101, 33.646335685700905, 51.75079495647569, 58.0715700634802, 40.08152068727785, 40.07834198735969, 32.178825492964826, 37.13988752821537, 61.14732373793652, 62.19361787616062, 18.87863176317971, 51.66116626954448, 7.716408393438067, 46.57703486158983, 27.624252564244067, 52.56896792685466, 96.73654629388271, 70.51917507675604, 70.88961398689557, 13.435363568775704, 29.896804775186325, 46.448419187991604, 52.558049147079146, 65.72262632840045, 80.7081879882318, 62.32994262739193, 81.7099331176993, 24.54347654297705, 85.40471286715805, 17.848433530308526, 2.7522650535389936, 46.0152862386255, 30.36799789298017, 62.445021937629264, 50.90361531230572, 88.19481599699215, 39.58233695139077, 55.91186271102641, 54.00540735731322, 41.914578699542204, 56.73721589054884, 25.71909476385006, 44.15733565175931, 34.316446686344506, 48.61317696858073, 66.17016426881884, 3.4275168979892356, 25.388526977134017, 55.31273818036356, 99.31859596820908, 86.40150322184053, 58.28762828406025, 48.25280288754749, 57.90542739177133, 92.48072349369144, 78.23826444505276, 76.68726852604853, 39.135743547284164, 38.52539095595869, 15.942342098931316, 92.44290784120868, 22.95127380173646, 26.625247171026835, 31.175782289561916, 83.51175760146154]]
        # ======================================================================
        # 11
        # omega = [[-0.1221193914430605, -0.25392161046154627, -0.4162959558123649, 0.396091630563459, -0.949346768366838, 0.7844483326946463, -0.9092041830777082, -0.7920669033615597, -0.9970719203770082, 0.9322830164995224, 0.8795625933329647, -0.5395740992568494, 0.12492901428942016, 0.695362957976315, 0.7048917912044079, -0.3725916861278422, -0.9465371496217887, 0.8874552805571518, 0.553525662742437, -0.18095477893366363, -0.3376983425085631, 0.6018074987250197, -0.12151775436996837, 0.05096381925197391, 0.948566440132955, 0.9274431128878105, -0.2936972192237748, -0.9431404058882129, 0.8960168487529394, -0.2620785953143854, -0.851978483011161, -0.4690129858795, -0.7544716012124935, 0.8373297921741905, -0.9910751085333855, 0.5709495705662404, -0.8656414786587467, 0.2324124012828115, 0.791496385587731, 0.8469347765529509, -0.672572049343247, -0.50037849953837, 0.19787706968671737, 0.5042207684319737, 0.9996470807832967, 0.8609075097751151, -0.7353316961182299, 0.8783396923097122, 0.014860384047297615, -0.7044563811473066, -0.5992629646821048, -0.8919836914616301, 0.3643288117460022, 0.9430169909643349, 0.5296796801609107, 0.24886305065365377, 0.3155932680149989, 0.992248245687698, -0.41708022468456174, -0.38953572196988584, -0.6142055411734646, 0.45614809509056475, 0.019702625599855805, -0.9782034930268051, 0.3132822984383399, -0.10563878869593668, -0.8048127212303136, 0.37412187712996237, -0.7661328401214658, -0.20379692803623684, -0.5663859465818626, -0.1429706248696765, 0.012714647785781796, 0.7933609602997238, -0.8756478494399706, -0.2998335032130608, 0.3514660312708908, -0.40225234618357875, 0.6931241121728839, -0.2449156895744744, 0.4603367257233453, 0.7487727950353329, -0.7760609176248638, -0.9840176924300383, 0.40486870627470406, 0.6376725576360753, -0.4957247782496803, 0.06819748596287156, -0.4342129407903188, 0.6732704938626959, 0.1444452519602537, -0.7847788115367398, 0.7725329231421545, -0.8440590564941424, 0.09429084522729814, -0.7455909561292127, -0.21203338480549538, -0.030514931031481396, 0.5839559118636612, 0.7598473469387017]]
        # mu = [[23.563837789306096, 52.45169508145242, 69.34662176249547, 23.91032666381725, 58.49371610465792, 68.14792757761234, 10.897161744252559, 89.64493381047872, 27.46612065560541, 7.5632790463146415, 3.3262381187451617, 54.434006636978395, 22.992115440028094, 55.446489187004424, 35.396169940811205, 66.85353971685501, 51.25254634202007, 68.61809982944278, 21.624312097763642, 82.1398341998496, 47.85403610546383, 40.585799459791424, 81.68897302381085, 83.53430361210887, 46.28292877619693, 44.58909714988347, 90.46677618823814, 57.25587847403245, 20.894230402828352, 31.285720038743936, 99.07796967725406, 91.28525456516674, 10.603944334998571, 78.6339953157281, 71.27994016108413, 48.02166394222888, 46.94591441129998, 45.47052602222143, 66.1115261662804, 92.50700817436993, 69.3927199085555, 28.622622106649644, 32.443160120308704, 67.17326038808399, 85.32099163699242, 29.472740127038087, 25.949792934972137, 27.97282826427162, 60.74001331989525, 28.19963922222113, 89.56674682536743, 84.49580292234378, 89.75024197375349, 26.96322900470215, 34.7685042843355, 62.73837662369893, 28.887889174636662, 99.17108414059152, 40.98288996677063, 96.04663660395694, 0.2797705982545251, 33.06686647431738, 42.66707891993426, 65.15968374406307, 35.1246165207016, 79.3897223061782, 27.538573640932917, 77.6728572610758, 97.97327792504801, 26.71740893118759, 73.73384367064452, 12.332564690383208, 29.963961586222744, 30.152332191402685, 30.33899000612149, 68.59247406530348, 12.335146598659774, 23.78608731377826, 16.209311260645364, 78.93680669314104, 83.45096711792242, 66.95305178147758, 99.45318542845754, 69.91991531620818, 1.0417975223098286, 10.709656881163044, 25.24124278278218, 25.249307748082483, 76.96134190253537, 26.43027021381942, 51.59238392229951, 76.62681290939514, 79.99950297337965, 24.767945243177103, 90.59671366821132, 74.46307555552141, 10.918708234505159, 40.1981781346052, 28.117426769106014, 90.8517412124946]]
        # ======================================================================
        #



        print(x)
        print(omega)
        print(mu)
        print("======================================================================")
        biotic = []
        perturb = []
        exp_temp = []
        add = 0



        system_state = np.zeros(SPECIES_K+ENV_VARS)

        Eg = ENV_START[0]

        # INIT ABUNDANCE
        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(Eg-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            system_state[s_i] = a_star
        #INIT ENV START
        for _ in range(ENV_VARS):
            system_state[SPECIES_K+_] = ENV_START[_]


        #INIT BIOTIC FORCE
        #biotic.append(0)
        #INIT PERTURBING FORCE
        perturb.append(0.2)
        #INIT EXPONENTIAL TEMP RUNOFF
        exp_temp.append(ENV_START[0])


        ALIVE_THRESHOLD=0
        results = [[] for _ in range(SPECIES_K+ENV_VARS)]
        times_steps=[]

        for step in np.arange(TIME_START, TIME_END, TIME_STEP):
            #print(step)
            times_steps.append(step)
            for _ in range(SPECIES_K+ENV_VARS):
                results[_].append(system_state[_])
            add = 0
            k1 = TIME_STEP * rates_of_change_system_state(system_state)
            k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
            k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
            add = 1
            k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)
            system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)
        ENV_VAR_ALIVE_ZERO_END = system_state[SPECIES_K+0]


        del perturb[-1]
        #del biotic[-1]
        del exp_temp[-1]



        fig, ax = plt.subplots()

        plt.title("Forces acting on the system", fontsize = 15)
        plt.xlabel('Time', fontsize = 15)
        plt.ylabel('Force', fontsize = 15)

        ax.plot(times_steps,perturb, 'b--',label = 'Perturbing Force')
        ax.plot(times_steps,biotic, 'g-',label = 'Biotic Force')

        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())


        ax.tick_params(which='both', width=1)
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=4)

        plt.xlim([0, 200])
        plt.ylim([-100, 50])

        plt.legend(prop={'size': 12}, loc='lower right')
        plt.tight_layout()
        plt.savefig(str(x)+".jpg")
        plt.show()



        fig, ax = plt.subplots()

        plt.title("System Temperature", fontsize = 15)
        plt.xlabel('Time', fontsize = 15)
        plt.ylabel('Temperature', fontsize = 15)

        ax.plot(times_steps,results[-1], 'r-',label = 'With alive species')
        ax.plot(times_steps,exp_temp, 'k--',label = 'Without alive species')

        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())


        ax.tick_params(which='both', width=1)
        ax.tick_params(which='major', length=7)
        ax.tick_params(which='minor', length=4)

        plt.xlim([0, 200])
        plt.ylim([0, 100])

        plt.legend(prop={'size': 12}, loc='lower right')
        plt.tight_layout()
        plt.savefig(str(x)+".jpg")
        plt.show()
