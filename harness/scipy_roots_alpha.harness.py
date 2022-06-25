import os
from multiprocessing import Process, Pool
import time
import numpy as np
import sys
from tqdm import tqdm

epoch_time = int(time.time())

def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    return(current_time)

K = 201
N = [5, 7, 10]

SAMPLE = 1

def run_it(Arg):
    Krun = Arg[0]
    Nrun = Arg[1]
    #print(Krun,Nrun)
    os.system("python3.10 " + os.getcwd() + "/experiments/scipy_roots_alpha.py " + str(Krun) + " " + str(Nrun) + " " + str(epoch_time))

if __name__ == '__main__':

    pool = Pool(processes=16)

    arguments_for_parallel = []

    for iter_k in range(1, K):
        for iter_n in N:
            for each_sample in range(SAMPLE):
                arguments_for_parallel.append((iter_k, iter_n))

    print("Starting: " + str(print_time()))

    results = []
    for result in tqdm(pool.imap_unordered(run_it, [_ for _ in arguments_for_parallel]), total=len(arguments_for_parallel)):
        results.append(result)

    print()
    print("Completed: "+str(print_time()))


