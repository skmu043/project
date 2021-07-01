import os
from multiprocessing import Process, Pool
import time
import numpy as np

epoch_time = int(time.time())

#print(epoch_time)'
def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(current_time)

SAMPLE = 100

def run_it(phi):
    #"e.g                                K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN = 2, OE = 5, LP_Z = (10 - 100), RUN_ID=epoch"
    print("Running with PHI : ", phi)
    os.system("python3.9 " + os.getcwd() + "/experiments/dyke_space.py 100 100 0 10 0 200 0.1 2 5 " + str(phi) + " " + str(epoch_time))

if __name__ == '__main__':

    pool = Pool(processes=8)

    for phi in np.arange(10, 91, 10):
        print_time()
        pool.map(run_it, [phi for x in range(SAMPLE)])
        print("Completed : " + str(phi))




