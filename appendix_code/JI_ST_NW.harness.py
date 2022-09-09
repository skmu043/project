import os
import time
import random
import shelve
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool


RUN_ID = int(time.time())                   # --- Unique Run ID
SAMPLE_SIZE = 100                           # --- Size of sample
SAMPLE_STEP = 1                             # --- Sample Step
SPECIES_K           = 100                   # --- Number of Species
RANGE_R             = 100                   # --- Essential Range
TIME_START          = 0                     # --- Start of Simulation
TIME_END            = 200                   # --- Length of Simulation
TIME_STEP           = 1                     # --- Time Steps
ENV_VARS            = 1                     # --- Number of Environment Variables
NICHE               = 5                     # --- Niche Size
SURVIVAL_THRESHOLD  = 0
ENV_START           = [50]
exp_name            = "JI_ST_NW.exp"

def init_shelve():
    data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
    shelve_file = data_directory + "/" + exp_name + ".data"
    os.mkdir(data_directory)
    s = shelve.open(shelve_file)

    try:

        s['RANGE_R']        = RANGE_R
        s['SPECIES_K']      = SPECIES_K
        s['TIME_START']     = TIME_START
        s['TIME_END']       = TIME_END
        s['TIME_STEP']      = TIME_STEP
        s['ENV_VARS']       = ENV_VARS
        s['exp_name']       = exp_name
        s['data_directory'] = data_directory
        s['shelve_file']    = shelve_file

    finally:
        s.close()

    return shelve_file

def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return(current_time)

def run_experiments(simulation_run_shelve):
    os.system("python3.10 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(simulation_run_shelve))

if __name__ == '__main__':

    shelve_files = []
    samples_sequence = np.arange(0, SAMPLE_SIZE, SAMPLE_STEP)
    print("Starting File Writes: " + print_time())

    for _ in tqdm(samples_sequence):

        omega               = [[random.uniform(-1, 1) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        mu                  = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]

        for start_temperature in np.arange (0,101, 5):
            for survival_threshold in np.arange (0,0.3, 0.2):
                for niche_size in [5, 10]:
                    if not (niche_size == 10 and survival_threshold == 0): # Math Domain Error : The Survival threshold cannot be zero
                        simulation_run_shelve = init_shelve()
                        shelve_files.append(simulation_run_shelve)
                        simulation_shelve = shelve.open(simulation_run_shelve)

                        try:
                            simulation_shelve['omega']              = omega
                            simulation_shelve['mu']                 = mu
                            simulation_shelve['ENV_START']          = start_temperature
                            simulation_shelve['SURVIVAL_THRESHOLD'] = float(str("{:.2f}".format(survival_threshold)))
                            simulation_shelve['NICHE']              = niche_size

                        finally:
                            simulation_shelve.close()

    print("Completed File Writes: " + print_time())

    print("Starting parallel Experiments: " + print_time())

    pool = Pool(processes=32)
    results = []

    for result in tqdm(pool.imap_unordered(run_experiments, [_ for _ in shelve_files]), total=len(shelve_files)):
        results.append(result)

    print("Completed parallel Experiments: " + print_time())
