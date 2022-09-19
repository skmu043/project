import sys, os
import shelve
from tqdm import tqdm

RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data_ST_NW_100'
data_archives = os.listdir(data_dr)

for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.niche.gaussian.exp.data")
    try:
        omega = s['omega']
        mu = s['mu']
        niche = s['NICHE']
        survival_threshold = s['SURVIVAL_THRESHOLD']

        data_input = (omega,mu,niche,survival_threshold)
        RESULT_DATA.append((data_input))

        if(data_input not in UNIQ_SAMPLES):
            UNIQ_SAMPLES.append(data_input)

            output_file = open("samples.file", "a")
            output_file.write(str(data_input))
            output_file.close()
    finally:
        s.close()


#=======================================================================================================================
def data_verification():

    for data_point in RESULT_DATA:
        if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
            UNIQ_SAMPLES.append((data_point[0],data_point[1]))


    DATA_VERIFICATION = []

    for uniq_ in UNIQ_SAMPLES:
        niche_stats = []
        survival_threshold_stats = []

        for data_point in RESULT_DATA:
            if uniq_ == ((data_point[0], data_point[1])):
                if([data_point[2],0] not in niche_stats):
                    niche_stats.append([data_point[2], 0])
                if([data_point[3],0] not in survival_threshold_stats):
                    survival_threshold_stats.append([data_point[3], 0])

        niche_stats.sort()
        survival_threshold_stats.sort()


        index_ = 0
        for niche_s in (niche_stats):
            for data_point in RESULT_DATA:
                if (uniq_ == ((data_point[0], data_point[1])) and data_point[2] == niche_s[0]):
                    niche_stats[index_][1] += 1

            index_ +=1

        index_ = 0
        for st_v in (survival_threshold_stats):
            for data_point in RESULT_DATA:
                if (uniq_ == ((data_point[0], data_point[1])) and data_point[3] == st_v[0]):
                    survival_threshold_stats[index_][1] += 1

            index_ +=1


        DATA_VERIFICATION.append([niche_stats,survival_threshold_stats])

    niche_check = []
    st_check = []

    print("DATA VALIDATION CHECK : ")
    print(all(i == DATA_VERIFICATION[0] for i in DATA_VERIFICATION))

    print("=====================")

data_verification()