import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor_core.data")

    try:

        RUN_ID = s['RUN_ID']

        number_alive_global_start = s['number_alive_global_start']
        number_alive_local_start = s['number_alive_local_start']
        number_alive_global_end = s['number_alive_global_end']
        number_alive_local_end = s['number_alive_local_end']

        number_alive_start = s['number_alive_start']
        number_alive_end = s['number_alive_end']

        # SUPER DATA STRUCTURE NEEDED

        #(Eg, El, number_alive_global_start, number_alive_local_start, number_alive_start, number_alive_global_end, number_alive_local_end, number_alive_end)

        #Eg = []
        #El = []
        #number_alive_global_start= []

        # Loops would be for item in Eg ...
        # plot Eg = x , El = y (the correct one) , heatmap += number_alive_global_start

        # Replicate for the ones below :

        #number_alive_local_start
        #number_alive_start
        #number_alive_global_end
        #number_alive_local_end
        #number_alive_end

    finally:
        s.close()

if __name__ == '__main__':

    print("Completed")
