import os
from multiprocessing import Process
import time
import numpy as np

epoch_time = int(time.time())

print(epoch_time)

def run_it(phi):
    print(phi)
    time.sleep(5)

if __name__ == '__main__':

    jobs=[]

    for phi in np.arange(0.1, 1, 0.1):
        print(phi)
        for _ in range(2):
            print(_, phi)
            jobs.append(Process(target=run_it,args=(phi,)))

    for job in jobs:
        job.start()

    for job in jobs:
        job.join()


#"Args: K, R, P, E, start, end, step, ENumber, Niche, PSR, Local Population Size"
#"e.g                                K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN = 2, OE = 5, LP_Z = (10 - 100)"
os.system("python3.9 " + os.getcwd() + "/experiments/dyke_space.py 100 100 0 10 0 200 0.1 2 5 10")
