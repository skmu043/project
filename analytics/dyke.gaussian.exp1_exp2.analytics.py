import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)
shel = shelve.open(data_dr + "/" + str(data_archives[0]) + "/dyke.gaussian.rk4.data")
RANGE_R = 0
times_step = 0

try:
    RANGE_R = shel['RANGE_R']
finally:
    shel.close()

plt.figure(figsize=(20,10))
plt.title('END ENV Var with 0 and 0.5 truncation for 50 start', fontsize=40)
plt.xlabel('Zero Truncation End Temperature', fontsize=20)
plt.ylabel('With Truncation End Temperature', fontsize=20)

plt.axhline(y=0, color='b', linestyle='-')
plt.axhline(y=100, color='b', linestyle='-')
plt.axvline(x=0, color='b', linestyle='-')
plt.axvline(x=100, color='b', linestyle='-')

plt.xlim([-50, 150])
plt.ylim([-50, 150])


for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.gaussian.rk4.data")

    try:

        ENV_VAR_ALIVE_ZERO_START = s['ENV_VAR_ALIVE_ZERO_START']
        ENV_VAR_ALIVE_ONE_START = s['ENV_VAR_ALIVE_ONE_START']
        ENV_VAR_ALIVE_ZERO_END = s['ENV_VAR_ALIVE_ZERO_END']
        ENV_VAR_ALIVE_ONE_END  = s['ENV_VAR_ALIVE_ONE_END']

        #print(ENV_VAR_ALIVE_ZERO_START, ENV_VAR_ALIVE_ZERO_END, "+ ", ENV_VAR_ALIVE_ONE_START, ENV_VAR_ALIVE_ONE_END)

        plt.scatter(ENV_VAR_ALIVE_ZERO_END,ENV_VAR_ALIVE_ONE_END)

    finally:
        s.close()

plt.show()