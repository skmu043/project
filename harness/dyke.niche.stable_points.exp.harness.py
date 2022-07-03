import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
import sys
from tqdm import tqdm

exp_name            = "dyke.niche.gaussian.exp.data"

RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)


for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.niche.gaussian.exp.data")
    try:
        omega = s['omega']
        mu = s['mu']
        RESULT_DATA.append((omega,mu))

    finally:
        s.close()


exp_name            = "dyke.niche.stable_points.exp"

for data_point in RESULT_DATA:
    if ((data_point[0],data_point[1])) not in UNIQ_SAMPLES:
        UNIQ_SAMPLES.append((data_point[0],data_point[1]))


def init_shelve():
    data_directory = str(os.getcwd())+"/data_stable/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
    shelve_file = data_directory + "/" + exp_name + ".data"

    os.mkdir(data_directory)
    s = shelve.open(shelve_file)

    try:
        s['exp_name'] = exp_name
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



    print("STARTING FILES: " + print_time())
    for each_sample in tqdm(UNIQ_SAMPLES):
        simulation_run_shelve = init_shelve()
        shelve_files.append(simulation_run_shelve)
        simulation_shelve = shelve.open(simulation_run_shelve)

        try:
            simulation_shelve['omega']              = each_sample[0]
            simulation_shelve['mu']                 = each_sample[1]

        finally:
            simulation_shelve.close()

    print("COMPLETED FILES: " + print_time())

    print("STARTING SIMULATIONS: " + print_time())
    pool = Pool(processes=16)

    results = []

    for result in tqdm(pool.imap_unordered(run_it, [_ for _ in shelve_files]), total=len(shelve_files)):
        results.append(result)

    #print(results)

    print("SIMULATIONS COMPLETED: " + print_time())

