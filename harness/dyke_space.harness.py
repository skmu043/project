import os
from multiprocessing import Process, Pool
import time
import numpy as np

epoch_time = int(time.time())

#print(epoch_time)

def run_it(one):
    #"Args: K, R, P, E, start, end, step, ENumber, Niche, PSR, Local Population Size"
    #"e.g                                K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN = 2, OE = 5, LP_Z = (10 - 100), RUN_ID=epoch"
    phi = 10
    os.system("python3.9 " + os.getcwd() + "/experiments/dyke_space.py 100 100 0 10 0 200 0.1 2 5 " + str(phi) + " " + str(epoch_time))

if __name__ == '__main__':

    jobs=[]

    pool = Pool(processes=4)
    #result = pool.apply_async(f, [10])
    #print result.get(timeout=1)
    pool.map(run_it, range(100))



    #for phi in np.arange(10, 100, 10):
        #print(phi)
        #for _ in range(100):
            #print(_, phi)

            #jobs.append(Process(target=run_it,args=(10,)))

    #for job in jobs:
        #job.start()

    #for job in jobs:
        #job.join()



