import random
import os
import shelve
from multiprocessing import Process, Pool
import time

# Generating ALL Parameters

SAMPLE_SIZE = int(10)
SAMPLE_STEP = int(50)
RUN_ID = int(time.time())

biotic_components_K = int(100)
essential_range_R = int(100)
external_perturbation_rate_P = int(0)
time_start = int(0)
time_end = int(200)
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
# End Generating Parameters

exp_name = "dyke_refactor"
data_directory = str(os.getcwd())+"/data/" + str(time.time()) + "." + str(random.randint(100, 999)) + "." + exp_name
shelve_file = data_directory + "/" + exp_name + ".data"
# Create Shelve to store parameters being sent to experiment run


def init_shelve():
    os.mkdir(data_directory)
    s = shelve.open(shelve_file)
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
        s['affects_w'] = affects_w
        s['optimum_condition_u'] = optimum_condition_u
        s['biotic_force_F'] = biotic_force_F

        s['exp_name'] = exp_name
        s['data_directory'] = data_directory
        s['shelve_file'] = shelve_file

    finally:
        s.close()


def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(current_time)


def run_it(phi):
    print("Running with PHI : ", phi)
    os.system("python3.9 " + os.getcwd() + "/experiments/" + exp_name + ".py " + str(RUN_ID))


if __name__ == '__main__':
    init_shelve()
    print_time()

    pool = Pool(processes=1)
    pool.map(run_it, [100 for x in range(SAMPLE_SIZE)])
