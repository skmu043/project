import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
import sys
from tqdm import tqdm

SAMPLE_SIZE = 50
SAMPLE_STEP = 1
RUN_ID = int(time.time())

SPECIES_K           = 100                   # ----------- Number of Biotic Components
RANGE_R             = 100                   # ----------- Essential Range
TIME_START          = 0                     # ----------- Start of Simulation
TIME_END            = 200                   # ----------- Length of Simulation
TIME_STEP           = 1                     # ----------- Time Step3
ENV_VARS            = 1                     # ----------- Number of Environment Variables
NICHE               = 5                     # ----------- Niche Size
LOCAL_SIZE          = 50                    # ----------- Local Population Size (%)
SURVIVAL_THRESHOLD  = 0
ENV_START           = [50]
exp_name            = "dyke.gaussian.exp1_exp2"

def init_shelve():
    data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
    shelve_file = data_directory + "/" + exp_name + ".data"

    os.mkdir(data_directory)
    s = shelve.open(shelve_file)
    #print(data_directory, shelve_file)

    try:

        s['RANGE_R']        = RANGE_R
        s['SPECIES_K']      = SPECIES_K
        s['TIME_START']     = TIME_START
        s['TIME_END']       = TIME_END
        s['TIME_STEP']      = TIME_STEP
        s['exp_name']       = exp_name
        s['data_directory'] = data_directory
        s['shelve_file']    = shelve_file
        #s['omega']          = omega
        #s['mu']             = mu
        s['ENV_VARS']       = ENV_VARS
        #s['NICHE'] = NICHE
        #s['SURVIVAL_THRESHOLD'] = SURVIVAL_THRESHOLD
        #s['ENV_START'] = ENV_START

    finally:
        s.close()

    return shelve_file

def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return(current_time)

#@jit
def run_it(simulation_run_shelve):

    # Main Experiment
    #print(simulation_run_shelve)
    os.system("python3.10 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(simulation_run_shelve))



if __name__ == '__main__':

    shelve_files = []

    samples_sequence = np.arange(0, SAMPLE_SIZE, SAMPLE_STEP)
    print("STARTING FILES: " + print_time())
    for _ in tqdm(samples_sequence):
        tqdm.write(str(_))
        omega               = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu                  = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

        for start_temperature in np.arange (5,100, 5):
            for survival_threshold in np.arange (0,1, 0.2):
                for niche_size in [5]:

                    #print(start_temperature, float(str("{:.2f}".format(survival_threshold))), niche_size)

                    simulation_run_shelve = init_shelve()
                    shelve_files.append(simulation_run_shelve)
                    #print("Creating: ",simulation_run_shelve)
                    simulation_shelve = shelve.open(simulation_run_shelve)

                    try:
                        simulation_shelve['omega']              = omega
                        simulation_shelve['mu']                 = mu
                        simulation_shelve['ENV_START']          = start_temperature
                        simulation_shelve['SURVIVAL_THRESHOLD'] = float(str("{:.2f}".format(survival_threshold)))
                        simulation_shelve['NICHE']              = niche_size

                    finally:
                        simulation_shelve.close()

    print("COMPLETED FILES: " + print_time())

    print("STARTING SIMULATIONS: " + print_time())
    pool = Pool(processes=16)

    #pool.map(run_it, [_ for _ in shelve_files])


    results = []

    for result in tqdm(pool.imap_unordered(run_it, [_ for _ in shelve_files]), total=len(shelve_files)):
        results.append(result)

    #print(results)

    print("SIMULATIONS COMPLETED: " + print_time())


#    data_dr = os.getcwd() + '/data'
#    data_archives = os.listdir(data_dr)

#    for file in data_archives:
#        sl = shelve.open(data_dr + "/" + str(file) + "/"+exp_name+".data")

#        try:

#            o = sl['omega']
#            m = sl['mu']
#            e = sl['ENV_START']
#            st = sl['SURVIVAL_THRESHOLD']
#            n = sl['NICHE']
#            print(o)
#            print(m)
#            print()

#        finally:
#            sl.close()
