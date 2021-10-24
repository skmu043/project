import os
from multiprocessing import Process, Pool
import time
import numpy as np
import sys

epoch_time = int(time.time())

def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(current_time)

K = 201
N = [5, 7, 10]

SAMPLE = 10

def run_it(K, N):
    Krun = K
    Nrun = N
    print(Krun,Nrun)
    os.system("python3.9 " + os.getcwd() + "/experiments/biotic_effect.py " + str(Krun) + " " + str(Nrun) + " " + str(epoch_time))

if __name__ == '__main__':


    os.system("python3.9 " + os.getcwd() + "/experiments/biotic_effect.py " + str(100) + " " + str(5) + " " + str(epoch_time))
    sys.exit()

    pool = Pool(processes=8)

    arguments_for_parallel = []

    for iter_k in range(1, K):
        for iter_n in N:
            for each_sample in range(SAMPLE):
                arguments_for_parallel.append((iter_k, iter_n))

    print(arguments_for_parallel)
    pool.starmap(run_it, arguments_for_parallel, 1)




