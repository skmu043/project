import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
import sys

# Generating ALL Parameters
SAMPLE_SIZE = 1
SAMPLE_STEP = 50
RUN_ID = int(time.time())

SPECIES_K   = 100                    # ----------- Number of Biotic Components
RANGE_R     = 100                   # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 0.1                   # ----------- Time Step
ENV_VARS    = 2                     # ----------- Number of Environment Variables
NICHE = 5                     # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0.5
ENV_START=[]
omega = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

local_population_index = []
uniq_k = []
for _ in range(int(LOCAL_SIZE/100 * SPECIES_K)):
    local_species = random.randint(0,SPECIES_K-1)
    while local_species in local_population_index:
        local_species = random.randint(0,SPECIES_K-1)
    local_population_index.append(local_species)
local_population_index.sort()


# Create Shelve to store parameters being sent to experiment run
exp_name = "dyke.refactor.rk4"

def init_shelve():
    data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
    shelve_file = data_directory + "/" + exp_name + ".data"

    os.mkdir(data_directory)
    s = shelve.open(shelve_file)
    print(data_directory, shelve_file)


    try:
        #s['SAMPLE_SIZE'] = SAMPLE_SIZE
        #s['SAMPLE_STEP'] = SAMPLE_STEP
        s['RUN_ID'] = RUN_ID

        s['SPECIES_K'] = SPECIES_K
        s['RANGE_R'] = RANGE_R
        s['TIME_START'] = TIME_START
        s['TIME_END'] = TIME_END
        s['TIME_STEP'] = TIME_STEP
        s['ENV_VARS'] = ENV_VARS
        s['NICHE'] = NICHE
        s['LOCAL_SIZE'] = LOCAL_SIZE
        s['ALIVE_THRESHOLD'] = ALIVE_THRESHOLD
        #s['ENV_START'] = ENV_START

        s['exp_name'] = exp_name
        s['data_directory'] = data_directory
        s['shelve_file'] = shelve_file

    finally:
        s.close()

    return shelve_file

def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(current_time)
print_time()


def run_it(simulation_run_shelve):

    # Main Experiment
    os.system("python3.9 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(simulation_run_shelve))


if __name__ == '__main__':

    shelve_files = []

    for Eg_temp in np.arange(1,RANGE_R,SAMPLE_STEP):
        for El_temp in np.arange(1,RANGE_R,SAMPLE_STEP):
            #print(Eg_temp, El_temp)

            simulation_run_shelve = init_shelve()
            shelve_files.append(simulation_run_shelve)
            print("InLoopCreates: ",simulation_run_shelve)
            simulation_shelve = shelve.open(simulation_run_shelve)

            ENV_START=[]
            ENV_START.append(Eg_temp)
            ENV_START.append(El_temp)

            try:
                simulation_shelve['omega'] = omega
                simulation_shelve['mu'] = mu
                simulation_shelve['local_population_index'] = local_population_index
                simulation_shelve['Eg'] = Eg_temp
                simulation_shelve['El'] = El_temp
                simulation_shelve['ENV_START'] = ENV_START

            finally:
                simulation_shelve.close()
            ENV_START.clear()

    print("===")
    #$for item in shelve_files:
     #   print(item)


    pool = Pool(processes=7)

    pool.map(run_it, [_ for _ in shelve_files])

    print("Completed")
