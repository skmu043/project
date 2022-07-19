import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
import sys
from tqdm import tqdm

exp_name            = "dyke.attractors"

data_dr = os.getcwd() + '/data_attractors'
data_archives = os.listdir(data_dr)

SAMPLE_SIZE = 7

def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return(current_time)

def run_it(simulation_run_shelve):
    os.system("python3.10 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(simulation_run_shelve))

if __name__ == '__main__':

    shelve_files = []

    print("STARTING FILES: " + print_time())

    exp_name_list = []

    for each_sample in range(SAMPLE_SIZE):
        exp_name_list.append(exp_name)

    print("COMPLETED FILES: " + print_time())

    print("STARTING SIMULATIONS: " + print_time())
    pool = Pool(processes=7)

    results = []

    for result in tqdm(pool.imap_unordered(run_it, [_ for _ in exp_name_list]), total=len(exp_name_list)):
        results.append(result)

    print("SIMULATIONS COMPLETED: " + print_time())

