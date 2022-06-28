import random
import os
import shelve
import time
from multiprocessing import Process, Pool
import numpy as np
import time
from matplotlib.gridspec import GridSpec

import sys
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

import statistics

SAMPLE_SIZE = 1
SAMPLE_STEP = 1
RUN_ID = int(time.time())

SPECIES_K   = 100                  # ----------- Number of Biotic Components
RANGE_R     = 100                  # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 1                   # ----------- Time Step3
ENV_VARS    = 1                     # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0
ENV_START=[90]

def abundance_loss(mu):

    truncation = 0.2
    step = 0.1

    ji_sum_total = 0
    ji_sum = [[] for _ in range(SPECIES_K)]
    for y in range(SPECIES_K):
        for x in np.arange (0, RANGE_R, step):
            ji_sum[y].append((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2)))))
        ji_sum_total+=sum(ji_sum[y])


    st_sum_total = 0
    st_sum = [[] for _ in range(SPECIES_K)]

    for y in range(SPECIES_K):
        for x in np.arange (0, RANGE_R, step):
            aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
            if(aliveness <= truncation):
                st_sum[y].append(0)
            else:
                st_sum[y].append(aliveness)
        st_sum_total += sum(st_sum[y])

    #print("Sum Abundance JI : ")
    #print(ji_sum_total)
    #print("Sum Abundance ST : ")
    #print(st_sum_total)

    return(ji_sum_total, st_sum_total)


if __name__ == '__main__':

    percent_lost = []

    for _ in range(100000):
        mu = [[random.uniform(0, RANGE_R) for _ in range(SPECIES_K)] for _ in range(ENV_VARS)]
        ji_st_ab = abundance_loss(mu)
        diff = ji_st_ab[0]-ji_st_ab[1]
        print(diff)
        loss = (diff)/ji_st_ab[0] * 100
        print(loss)
        percent_lost.append(loss)
        print("mean % loss")
        print(statistics.mean(percent_lost))
