import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
import sys

# Generating ALL Parameters
SAMPLE_SIZE = int(1)
SAMPLE_STEP = int(1)
RUN_ID = int(time.time())

biotic_components_K = int(50)
essential_range_R = int(100)
external_perturbation_rate_P = int(0)
time_start = int(0)
time_end = int(200)
time_step = float(1)
environment_components_N = int(2)
truncated_gaussian_ROUND = int(1)
niche_width = int(5)
local_population_size = int(100)
biotic_force_F = [0 for _ in range(environment_components_N)]

affects_w = [[] for _ in range(environment_components_N)]
for wi in range(environment_components_N):
    affects_w[wi] = [random.uniform(-1, 1) for _ in range(biotic_components_K)]
for wi in range(environment_components_N):
    while int(len(affects_w[wi])) != int(len(set(affects_w[wi]))):
        print("Duplicate w's detected: Regenerating ...")
        affects_w[wi].clear()
        affects_w[wi] = [random.uniform(-1, 1) for _ in range(biotic_components_K)]

optimum_condition_u = [[] for _ in range(environment_components_N)]
for ui in range(environment_components_N):
    ######################################################### BIG CHANGE HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #optimum_condition_u[ui] = [random.uniform(0, essential_range_R-1) for _ in range(biotic_components_K)]
    optimum_condition_u[ui] = [random.uniform(0, essential_range_R) for _ in range(biotic_components_K)]
for ui in range(environment_components_N):
    while int(len(optimum_condition_u[ui])) != int(len(set(optimum_condition_u[ui]))):
        print("Duplicate u's detected: Regenerating ...")
        optimum_condition_u[ui].clear()
        optimum_condition_u[ui] = [random.uniform(0, essential_range_R) for _ in range(biotic_components_K)]

for item in optimum_condition_u:
    print(item)

local_population_index = []
uniq_k = []
for x in range(int(local_population_size/100 * biotic_components_K)):
    one = random.randint(0,biotic_components_K-1)
    while one in uniq_k:
        one = random.randint(0,biotic_components_K-1)
    uniq_k.append(one)
    local_population_index.append(one)

local_population_index.sort()

# Create Shelve to store parameters being sent to experiment run
exp_name = "truncated_gaussian"

def init_shelve():
    data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
    shelve_file = data_directory + "/" + exp_name + ".data"

    os.mkdir(data_directory)
    s = shelve.open(shelve_file)
    print(data_directory, shelve_file)

    try:
        s['SAMPLE_SIZE'] = SAMPLE_SIZE
        s['SAMPLE_STEP'] = SAMPLE_STEP
        s['RUN_ID'] = RUN_ID
        s['biotic_components_K'] = biotic_components_K
        s['essential_range_R'] = essential_range_R
        s['external_perturbation_rate_P'] = external_perturbation_rate_P
        s['time_start'] = time_start
        s['time_end'] = time_end
        s['time_step'] = time_step
        s['environment_components_N'] = environment_components_N
        s['truncated_gaussian_ROUND'] = truncated_gaussian_ROUND
        s['niche_width'] = niche_width
        s['local_population_size'] = local_population_size
        s['biotic_force_F'] = biotic_force_F
        s['exp_name'] = exp_name
        s['data_directory'] = data_directory
        s['shelve_file'] = shelve_file

    finally:
        s.close()

    return shelve_file

def run_once(simulation_run_shelve):

    os.system("python3.9 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(simulation_run_shelve))

if __name__ == '__main__':

    shelve_files = []
    simulation_run_shelve = init_shelve()
    simulation_shelve = shelve.open(simulation_run_shelve)

    try:
        simulation_shelve['affects_w'] = affects_w
        simulation_shelve['optimum_condition_u'] = optimum_condition_u
        simulation_shelve['local_population_index'] = local_population_index
    finally:
        simulation_shelve.close()

    run_once(simulation_run_shelve)
