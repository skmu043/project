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

K       = 100        #Number of Biotic Components
R       = 100        #Essential Range (defines where Biotic Components can be present)
OE      = []         #Niche
start   = 0          #Time Start
end     = 300        #Time End
step    = 0.1        #Time Step
w       = []         #Affects Parameter (ranges between -1 and 1 for each K)
u       = []         #Ideal Temperature for species (between 0 and R -> the essential range)

# The following will change from being single to multiple to support multiple E
N       = 2          #Number of Environment Variables
Ei      = -20

E       = [random.uniform(-1,10) for _ in range(N)]
#E       = [50 for _ in range(N)]

#same start off for E
#E       = [Ei for _ in range(N)]

Pi      = 0
P       = [Pi for _ in range(N)]
F       = [Pi for _ in range(N)]
#individual P rates
Px      = [random.uniform(0.1,0.3) for _ in range(N)]
#zero P rates
#Px      = [0 for _ in range(N)]

Et      = Ei         #Temperature without Biotic Force

#OE = [random.uniform(3,10) for _ in range(K)] #Switches between same sized Niches to different sized ones
OE = [5 for _ in range(K)]

#populates affects values
w = [[] for _ in range(N)]
for wi in range(N):
    w[wi] = [random.uniform(-1,1) for _ in range(K)]

#print(w)
for wi in range(N):
    while (int(len(w[wi])) != int(len(set(w[wi])))):
        print("Duplicate w's detected: Regenerating ...")
        w[wi].clear()
        w[wi] = [random.uniform(-1,1) for _ in range(K)]

#print(w)

#w = [0.5342271292313454, -0.20219795629536352, -0.44480834167189354, 0.22569487814594358, 0.746963728629346, 0.6175325100969504, -0.7945015070961714, -0.34018517438798135, -0.30369966581074515, 0.7910181158182961, -0.8785502663905314, 0.1303023728474253, -0.46157360594593677, -0.6730578750616567, -0.5320788094747058, 0.3376222547693799, -0.566209933186111, 0.1323326828728253, 0.671260375790319, 0.842586189610675, 0.22088077662193362, 0.6938115823367961, -0.26650298525774563, -0.9287537472184237, 0.008187773679745769, -0.32827674194154355, -0.4669427128677155, 0.33626407608453146, 0.021921885191599433, 0.5307543086114925, -0.5457348859903011, 0.9051977392015649, -0.1812621059222359, 0.29479991810834183, -0.43832179499171464, -0.5727671925482505, -0.22470436322927534, -0.8825523065836953, 0.6993841890376307, 0.012855365095307159, -0.7206439819061259, -0.32547564179472266, -0.15061958357567873, 0.38797658079751307, -0.9233198482143754, -0.03491465819950679, 0.4906750396170554, -0.1971338922348742, -0.49200574618643356, 0.6765869421626849, 0.7001386902159443, -0.05496321044367192, -0.2736181880924944, -0.7801912949328946, -0.5996113503206131, -0.7859665002768352, -0.037655705240655646, -0.4940315977094547, -0.16212529971201062, 0.3706800324065487, -0.5442133942209184, -0.594034431637011, -0.29814513702218504, -0.49590937738482954, 0.9111264058260566, 0.3738124323744465, 0.7887849059439367, 0.9691569067086034, 0.9925775673328701, 0.7817628857069394, -0.13111752712541702, -0.9275418680823477, 0.013428250501857253, 0.7559737226311734, 0.6848054000476091, -0.14745807424988877, -0.29536914081469523, -0.6531323767279722, -0.6438351088045227, 0.1844680697669785, -0.376087115503813, 0.03967810787769488, -0.09549136842652017, 0.2791317531986153, -0.6595896382600808, -0.37739341864542997, 0.4835077411924358, 0.047304126640932864, 0.9488183568883797, -0.13134627165778934, 0.7581147679818867, -0.7949874553975236, -0.7320962687423667, 0.5371391158099199, -0.39053599834586117, -0.5428923946311217, -0.24624524264681935, 0.1423221391012408, -0.0899201441097941, 0.7915526016336334, -0.3496883694859061, 0.9519627788713105, 0.9753484727017931, 0.9166501702220766, -0.9031056698957745, -0.3113567321206785, 0.4560908739109024, 0.7336780404219241, 0.5255736924910277, 0.27464586014224834, 0.5708863936976487, -0.7156787868190175, 0.5567899085567991, 0.5774069428782673, -0.9504969012398092, 0.40388924236290213, -0.10061822441554913, 0.7538847194822476, 0.8950600330284746, -0.6046144967815015, -0.5046375051253382, 0.7216102984844905, -0.9053595194327584, 0.9304655588549418, 0.12931934946592283, -0.9223192014101183, -0.1626171882874341, 0.5342742955964839, 0.11392046519240151, 0.9157226619554297, 0.2601726967465383, 0.07724360525606278, -0.4133677095472652, 0.306043230178606, -0.15223615452138284, -0.3451544020350483, 0.39386446065772085, 0.435558476624174, -0.15104641708964972, -0.42012852888374885, -0.5035216227678259, -0.47614461569111644, 0.6463416634774639, 0.88302808108027, -0.5026724673780263, -0.5187768733156659, 0.14925561746703164, 0.4096739023942104, -0.7622148258032011, -0.22982304276750676, -0.7961716151380624, -0.08385564667741607, 0.5512336696660203, 0.40959836128315197, 0.8042297204681683, 0.821514821198801, -0.609178845039215, -0.7782385606171622, -0.733763785679733, 0.5432801756195784, -0.22944625216548786, 0.5076874060494236, 0.5183575876876672, 0.3459492337023422, -0.852859975816161, -0.5366670279096966, 0.17329118019196188, 0.7667611123026514, 0.00017026855997315238, 0.10431007370377965, -0.4044283489994105, 0.42978886148536777, -0.26955695206121, 0.7511147071983308, -0.8138965160701086, -0.10217816022068482, 0.3634901941184365, -0.28073803178566026, -0.9779321765430593, 0.9600278589280213, 0.9393715629166013, -0.7577028166922419, 0.5539023534881347, -0.31937213691148947, 0.0032107195959236723, -0.03534970363279233, -0.7426344035319616, 0.24678790810569629, -0.18165412113421042, 0.7048083099539366, 0.8393262143111286, -0.36265401005659426, 0.10026906923927936, 0.6450186663649538, -0.048155215286273645, 0.9000322412479285, -0.14293308717437547, 0.6396897844975524, 0.6446977622979264, -0.7998500298440485, -0.004580647767524315, -0.32710886270318773, 0.6336278388465235, 0.8383682799848655, 0.9268526967376276, -0.10720626573508962, -0.5423269125095922, -0.7966456583689183, -0.8905716724673609, -0.9256489170990736, -0.9224867175043097, 0.25284514268525826, 0.7881521666923503, -0.8769319957961677, 0.5804525141184766, -0.4368194289524123, 0.6673031479114742, -0.35612951584249397, -0.011800786345770309, -0.259701175326984, 0.6121758548002372, 0.067722710826436, -0.3172343719580284, -0.5839745355481862, -0.17153051693644317, -0.5358956137528994, 0.8287157908768688, 0.42777895638756713, -0.07288480295219535, 0.5257678500312841, 0.861581282137261, 0.03558013040594665, 0.818537371056038, 0.9824112262720719, -0.7751772644078729, -0.20544597384607477, -0.12162171074913686, 0.663229376786457, 0.22788323245830533, 0.5997513844950413, 0.21980115788231958, 0.8747143263003339, 0.7161672695157137, 0.5760036671802697, -0.9614163375242208, -0.9730530095271079, 0.8205108159335142, 0.9401878569099271, 0.07902919762900606, 0.016190418247956417, 0.07900784426435803, -0.5510611093623599, -0.06300174171181028, 0.07442816381478479, -0.4497094591526194, -0.06788415462588038, -0.9804999946440824, 0.9480902456957039, -0.15295650613734524, 0.09597599931708567, -0.6367425746408071, 0.8414605305082392, 0.7793278953211908, -0.9040809638782408, 0.49856431148370683, 0.291045122882182, -0.7618070857165566, -0.3431279511558287, -0.20941023118466395, 0.17915166040727826, 0.5522404180837075, -0.49201722974861073, -0.8726195748448864, -0.43062176153265574, -0.28613452641990644, 0.9499794657558025, -0.7982475905792359, 0.35155823542216225, 0.4307636927570513, -0.5058277012493251, -0.40623289275259866, 0.5941645418511281, 0.9198648995852241, 0.8920934744827449, 0.3687076524976083, 0.986988503085761, -0.4475520872080152, -0.6619648413069221, -0.26371084558137303, -0.19920429714744348, -0.7663222899753674, 0.5610490903746146, 0.8708779974452556, -0.1639301272599505, 0.717445713423779, -0.9799323173389047, 0.9916149702583685, 0.9274354808361243, 0.32016364081766535, -0.8466808256063403]
# #populates optimum growing temperatures [e.g 78.213423]
u = [random.uniform(0, R) for _ in range(K)]

while(int(len(u)) != int(len(set(u)))):
    print("Duplicate u's detected: Regenerating ...")
    u.clear()
    u = [random.uniform(0, R) for _ in range(K)]

#print(u)

#u = [1.9450593597140076, 86.35718651062018, 86.72638325100026, 12.729322580032, 57.16938733306887, 46.09789040038521, 98.63971722872468, 90.04427417879262, 78.62023629634596, 65.8962035318393, 55.85499976874695, 9.487715695190712, 2.3443574584210425, 45.71662599462844, 84.8839215011532, 2.8451584634601734, 84.21130987885158, 93.94615303478983, 95.37363778400318, 63.904021203450235, 49.72199553534442, 40.11224992811777, 4.648715035112472, 18.24048527090991, 11.872311160254444, 44.02712903700103, 10.642631784289058, 9.714403865223197, 56.86153216490566, 25.039863106686667, 12.479711612152688, 26.542718248540165, 12.046183150972622, 87.80787901132331, 56.95953061144038, 65.05207456401813, 12.490117909715314, 37.809283809980855, 37.90978137287804, 49.796627133260266, 83.18092392473521, 65.83599387199264, 55.366484549558656, 67.13956218510717, 45.75388872514773, 36.971319039811455, 60.42873207148622, 43.732521718993354, 2.303371173469415, 35.841954514917205, 3.039399306907453, 26.201883032962957, 53.42046238629629, 34.053767651158374, 52.77636429374798, 72.23802635656281, 70.51508379012483, 58.7898893108219, 42.88101578763251, 4.201085828747009, 13.167650529870055, 68.67203522374702, 47.79137311881421, 96.35780426172501, 80.30599018611181, 80.53216745728099, 48.4566659139069, 64.41107021625699, 73.22908999624744, 63.78515752541387, 68.019953942972, 46.95330029987771, 44.54246326249922, 58.029454883459906, 49.24367973816381, 57.93625253714198, 20.10953275651277, 9.894593420429798, 37.82794789521729, 0.322695994484512, 19.430916067054994, 90.79049685068846, 49.20598902680604, 63.06163076214929, 19.671055110153045, 40.323052754292256, 1.272717426670722, 20.464526275161777, 48.319959842316564, 54.242074112001085, 77.36033680183787, 60.23801215137363, 61.24507324294958, 63.39377412716692, 91.75789946422248, 55.52507145374377, 79.8162920335649, 32.73316669367009, 73.2315063470369, 40.60506680022681, 15.243052759952514, 40.41983076817905, 85.65389619459532, 41.646128265210876, 95.08828403635526, 87.03422753711787, 65.86213574030539, 2.660562738965766, 65.92028458103928, 46.019017482072485, 10.859251939561343, 56.935229274800946, 74.94659035503312, 64.69475007979803, 72.64693408355673, 4.558598812936909, 36.99760614642863, 43.056269902334485, 20.51775683063567, 40.77790208125593, 46.11308057819142, 64.50402716622892, 60.20513542368875, 70.48888252501439, 23.79642884634363, 87.67905598948677, 56.53271648460052, 21.350204934463367, 30.38959955731978, 32.75044027494071, 79.00150724522202, 38.221330005036904, 64.67888022009377, 41.79995027879689, 25.756888415913636, 84.84158263983018, 68.0496487261751, 16.014016141152776, 68.9518579359388, 66.95588231344915, 53.61630017080054, 60.242645182137345, 99.30922138353559, 14.326747680534325, 25.959787197538287, 48.75826639546565, 38.0574051158226, 9.930276980232245, 0.3206207062279032, 33.68457374879843, 91.38126378217068, 6.760555911793464, 32.78950554508609, 44.35589419805366, 84.73253854489047, 0.02498255186065279, 47.70270785610856, 91.11053898307166, 1.6210835826729997, 87.21586080111942, 59.05068075434225, 17.98207083156731, 90.76381896825644, 50.3789655255068, 9.812410479435618, 48.33203513581768, 90.04292933627104, 26.924440439026686, 32.3894236632592, 24.376126075965676, 70.39037483780074, 71.62107984046486, 59.06979313292674, 56.3461839323006, 86.34466646043045, 80.6169397595535, 83.11443025820671, 69.18507338653444, 38.41738395592643, 65.17502404211956, 54.874089819465624, 76.81711358181595, 63.4788806206172, 29.61495628507477, 14.612653930750508, 5.047365099193646, 14.080968257330639, 42.17668816752729, 75.14465089342237, 55.84340261621533, 26.048314616523626, 85.13203132940637, 23.01260813335384, 38.90524909223549, 10.943890258166466, 47.787323608262014, 94.52590328617532, 14.180763451150401, 48.71601589018552, 39.614164518184936, 60.031628701643015, 65.70719341725308, 74.04897556998726, 32.72914347879894, 93.35387358970014, 93.79927098582985, 19.5438971215253, 89.48103779527085, 72.38129806107399, 80.01534752994786, 13.977909353774276, 99.67030499680841, 11.839658637897088, 66.15000717948308, 29.650057479008208, 77.14482221516327, 98.08585514617435, 75.52365904994492, 33.41073782695999, 91.7695933794295, 11.702679741866007, 82.74847717324974, 95.19389110535464, 35.02868474624117, 5.649952906829913, 70.73674475673307, 36.45968304006374, 93.44229362524901, 87.7040289907413, 95.94324895529603, 77.27458798141112, 48.86645339837975, 57.49445518370179, 1.9008632804904058, 10.18321700173489, 12.917120671716143, 39.46344074645612, 72.83824010591961, 16.712675768382468, 71.37337854607553, 49.63430571945982, 29.81869518968283, 14.331063900336648, 21.197614120807273, 89.72061268174404, 37.67072725798572, 73.16535536430908, 98.3338797586508, 53.898289916782005, 36.87340465027905, 38.84719417368838, 86.77579018302232, 90.99206235316028, 28.352664657560045, 97.44361637580732, 53.8371257508019, 99.73090496689811, 83.82865741596738, 13.779553809576795, 72.54694917437176, 52.425214205594315, 87.879583230042, 37.61715326965109, 97.8130758558271, 32.72911819520801, 62.10951150989129, 56.00578891077931, 51.52722356398869, 64.79292454014788, 3.4612043554667515, 35.47053425349178, 42.08795503408712, 91.56357518924526, 44.656832250279436, 56.91486802968294, 38.81533510583752, 74.6043942399619, 97.2480658926733, 74.27838439726597, 31.959847330204237, 27.771916701635913, 40.48653119940704, 75.70828095428831, 76.45273830657942, 4.177600535343872, 77.69543424498691, 54.08497108161189, 31.328332564078785, 95.53249866798689, 91.39782885570665, 64.34513166063243, 24.390275903087165, 16.099526004745435, 9.61531579870256, 13.06497832442879, 99.69023769566138, 45.03505313512417, 1.8052961059840555, 14.277742184905817, 47.208112156764074]

#alpha = [[] for _ in range(N)]
#for ai in range(N):
#    alpha[ai] = [[] for _ in range(K)] #abundance value for a species

alpha = [[] for _ in range(K)]

# Future Fix - this calculates abundance for only the first Environment Variable

for _ in range(K):
    al = []
    for ai in range(N):
        al.append((math.e) ** ((-1) * (((abs((E[ai])-u[_])) ** 2) / (2*(OE[_]**2)))))
    alpha[_].append(np.prod(al))

rF = [[] for _ in range(N)]         #Biotic Force Values
rP = [[] for _ in range(N)]         #Perturbation Values
rE = [[] for _ in range(N)]         #Temperature with Biotic Force Values
rEt = []                            #Temperature without Biotic Force Values

#Abundance values over time
rAx = [[] for x in range(K)]

#Abundance values over time scaled up by R (Essential Range)
rAxR = [[] for x in range(K)]
#Tracks time (steps accumulation)
time = []

biotic_force = [[] for _ in range(N)]
for i in range(N):
    biotic_force[i] = [[] for _ in range(K)]
temperatures = []

#plot abundance of species over temperature
def plot_alphas():

    for x in np.arange (0, R, step):
        temperatures.append(x)

    for i in range(N):
        for y in range(K):
            for x in np.arange (0, R, step):
                biotic_force[i][y].append((math.e) ** ((-1) * (((abs(x-u[y])) ** 2) / (2*(OE[y]**2)))) * w[i][y])

    plt.figure(figsize=(25,15))
    plt.title('Biotic Force over Time', fontsize=40)
    plt.xlabel('Temperature', fontsize=40)
    plt.ylabel('biotic force (a * w)', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    for i in range(N):
        for _ in range(K):
            plt.plot(temperatures,biotic_force[i][_])
    for i in range(N):
        plt.plot(temperatures,np.sum((np.array(biotic_force[i], dtype=float)), axis=0), lw=4)

    plt.show()

local_population_size = 20
local_population_index = []
for x in range(local_population_size):
    local_population_index.append(random.randint(0,K-1))

print(local_population_index)

def update(step):
    global F, P, E, Et, rF, rP, rE, rEt, u, w

    fSUM = [0 for _ in range(N)]
    alpha_time_scale = 1
    temperature_time_scale = 0.2

    for _ in range(K):
        al = []
        for ei in range(N):
            #print(_)
            al.append((math.e) ** ((-1) * (((abs((E[ei])-u[_])) ** 2) / (2*(OE[_]**2)))))
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

    #fSUM = fSUM + (alpha[_][-1] * w[_])

    #abundance directly on the graph
    #alpha[_] = al
    #rAx[_].append(alpha[_])
    #fSUM = fSUM + (alpha[_] * w[_]) # Fixed


    #rAx[_].append(alpha[_][-1])
    #rAxR[_].append(alpha[_][-1] * R)
    for _ in range(K):
        for ei in range(N):
            fSUM[ei] = fSUM[ei] + (alpha[_][-1] * w[ei][_])

    #abundance directly on the graph
    #alpha[_] = al
    #rAx[_].append(alpha[_])
    #fSUM = fSUM + (alpha[_] * w[_]) # Fixed

    for ei in range(N):
        F[ei] = fSUM[ei] * 10
        #F = F + (fSUM * step)

        P[ei] = P[ei] + (Px[ei] * step)

        #P = 0
        #F = fSUM                  [Explore the linear increase for P]
        #P = P + (step/3.5)
        newE = E[ei] + (((P[ei] + F[ei]) * step))
        # E is old E and newE has the current value
        E[ei] = E[ei] + ((newE-E[ei]) * temperature_time_scale)

        # E not P ! This is the Temperature !
        # Incorrect one Et = Et + P
        # F becomes 0 - no biotic force as no biota
        Et = Et + ((P[ei] + 0) * step)

        rF[ei].append(F[ei])
        rP[ei].append(P[ei])
        rE[ei].append(E[ei])
    rEt.append(Et)

#plot affects values for each species
def plot_w():

    plt.figure(figsize=(20,5))
    plt.ylim(-1, 1)
    plt.title('Affects values for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Affects', fontsize=18)
    for i in range(N):
        plt.plot(w[i], 'k.', label='w')
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
        plt.plot(time,rAxR[x])

    plt.plot(time,rE[0], label='E Global', linewidth=6)
    plt.plot(time,rE[1], label='E Local', linewidth=6)

    plt.legend(loc=5, prop={'size': 30})
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
        if(w[i][x]==0):
            plt.plot(time,rAx[x],'k-')
        if(w[i][x]<0):
            plt.plot(time,rAx[x],'b-')
        if(w[i][x]>0):
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
    for i in range(N):
        plt.plot(time,rF[i], 'g-', label='F')
        plt.plot(time,rP[i], 'b--', label='P')
    plt.show()

#plot temperature value over time
def plot_e():
    plt.figure(figsize=(20,10))
    plt.title('Environment Variable E', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for i in range(N):
        plt.plot(time,rE[i], 'r-', label='E')
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
    for i in range(N):
        plt.plot(time,rF[i], 'g-', label = 'biotic force')
        plt.plot(time,rP[i], 'b--', label = 'perturbing force(rate)')
        plt.plot(time,rE[i], 'r-',label = 'temperature')
    plt.plot(time,rEt, 'k.',label = 'temperature (without biota)')
    #plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)
    plt.show()


def plot_ep(LP):
    plt.figure(figsize=(20,10))
    plt.title('Simulation Values over Time', fontsize=40)
    plt.xlabel('Time', fontsize=40)
    plt.ylabel('Environment Variables E', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(0, R)
    plt.xlim(0, end)
    for i in range(N):
        plt.plot(time,rE[i], 'r--',label = 'temperature')
    plt.axvline(x=LP)
    plt.show()

sys.stdout.write("[%s]" % (" " * K))
sys.stdout.flush()
sys.stdout.write("\b" * (K+1))


if __name__ == '__main__':

    LP = 50

    for xtime in np.arange (start, end, step):
        update(step)
        time.append(xtime)
        #print(xtime)

        if(xtime % 1 == 0):
            sys.stdout.write("-")
            sys.stdout.flush()
        #f(xtime == LP):
        #   for x in range(N):
        #       P[x] -= 15

        #if(xtime == 100):
        #    for x in range(N):
        #        P[x] -= 10
    

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
    #plot_ep(LP)               #plot temperature and large P
