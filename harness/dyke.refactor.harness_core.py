import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np

# Generating ALL Parameters
SAMPLE_SIZE = int(1)
SAMPLE_STEP = int(1)
RUN_ID = int(time.time())

biotic_components_K = int(100)
essential_range_R = int(100)
external_perturbation_rate_P = int(0)
time_start = int(0)
time_end = int(100)
time_step = float(0.1)
environment_components_N = int(2)
truncated_gaussian_ROUND = int(10)
niche_width = int(5)
local_population_size = int(10)
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
    optimum_condition_u[ui] = [random.uniform(0, essential_range_R) for _ in range(biotic_components_K)]
for ui in range(environment_components_N):
    while int(len(optimum_condition_u[ui])) != int(len(set(optimum_condition_u[ui]))):
        print("Duplicate u's detected: Regenerating ...")
        optimum_condition_u[ui].clear()
        optimum_condition_u[ui] = [random.uniform(0, essential_range_R) for _ in range(biotic_components_K)]

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
exp_name = "dyke.refactor_core"
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
shelve_file = data_directory + "/" + exp_name + ".data"

# Store source shelve file to read common
control_data = str(os.getcwd())+"/data/" + exp_name + "_control"
control_shelve_file = control_data + "/" + exp_name + ".control.data"

print("Control Shelve: "+control_shelve_file)

if(not os.path.isdir(control_data)):
    print("gate one")
    if(not os.path.exists(control_shelve_file)):
        os.mkdir(control_data)
        control_shelve = shelve.open(control_shelve_file)
        print("WRITING !!!")
        try:
            control_shelve['affects_w'] = affects_w
            control_shelve['optimum_condition_u'] = optimum_condition_u
            control_shelve['local_population_index'] = local_population_index
        finally:
            control_shelve.close()
            print("CLOSING SHELVE")

print("Control Shelve End")
def init_shelve():

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
        #s['affects_w'] = affects_w
        #s['optimum_condition_u'] = optimum_condition_u
        s['biotic_force_F'] = biotic_force_F
        #s['local_population_index'] = local_population_index

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


def run_it(global_local_start_temp):
    #open control shelve and pack common vars for individual simulation runs
    control_shelve = shelve.open(control_shelve_file)

    simulation_run_shelve = init_shelve()
    simulation_shelve = shelve.open(simulation_run_shelve)


    try:
        simulation_shelve['affects_w'] = control_shelve['affects_w']
        simulation_shelve['optimum_condition_u'] = control_shelve['optimum_condition_u']
        simulation_shelve['local_population_index'] = control_shelve['local_population_index']
        simulation_shelve['global_start_temp'] = global_local_start_temp[0]
        simulation_shelve['local_start_temp'] = global_local_start_temp[1]
    finally:
        simulation_shelve.close()
        control_shelve.close()

    os.system("python3.9 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(simulation_run_shelve))

if __name__ == '__main__':



    simulation_start_temperatures = []

    for Eg_temp in np.arange(1,essential_range_R,SAMPLE_STEP):
        for El_temp in np.arange(1,essential_range_R,SAMPLE_STEP):
            simulation_start_temperatures.append((Eg_temp, El_temp))

    print(simulation_start_temperatures)

    pool = Pool(processes=1)
    pool.map(run_it, [x for x in simulation_start_temperatures])

    #     for Eg_temp in np.arange(1,100,SAMPLE_STEP):
    #         for El_temp in np.arange(1,100,SAMPLE_STEP):
    #             print("Init : ", Eg_temp, El_temp)
    #             > Open Shelve and update and Close Shelve note run_it
    #             is its own run so open run_it's local shelve and update that
    # which creates the problem of run_it using different w and u values
    # so solution could be for this one shelve, generate the start pairs and then iterate ...
    # start_trajectories = ((1,1, shelve_file), (5,5, shelve_file),
    # (10,10, shelve_file), (15.15, shelve_file), (20,20, shelve_file))
    # pool.map (run_it, start_trajectories)
    # run_it(start_traj)
    # 
