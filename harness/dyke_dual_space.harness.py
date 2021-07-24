import os
from multiprocessing import Process, Pool
import time
import numpy as np

# ((L1)(L2)) -> entire thing is affected by Eg
# (L1) affected by EL1 * Eg
# (L2) affected by EL2 * Eg


epoch_time = int(time.time())

#print(epoch_time)'
def print_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S", t)
    print(current_time)

SAMPLE = 100

K = 100

def run_it(LP_Zx):
    #"e.g                                K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN = 2, OE = 5, LP_Z = (10 - 100), RUN_ID=epoch"
    LP1_size = LP_Zx[0]
    LP2_size = LP_Zx[1]

    print("Running : ", LP_Zx)

    os.system("python3.9 " + os.getcwd() + "/experiments/dyke_dual_space.py "+ str(K) +" 100 0 10 0 200 0.1 3 5 " + str(LP1_size) + " " + str(LP2_size) + " " + str(epoch_time))

if __name__ == '__main__':

    pool = Pool(processes=8)

    LP_Zx = []

    for perc in np.arange(10, 51, 10):
        LP1_size = int(perc/100 * K)
        LP2_size = int(K - LP1_size)
        for sn in range(SAMPLE):
            LP_Zx.append((LP1_size, LP2_size))
            #print_time()
            #print(LP1_size, LP2_size)

    pool.map(run_it, LP_Zx)



