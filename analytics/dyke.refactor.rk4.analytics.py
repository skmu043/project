import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)
shel = shelve.open(data_dr + "/" + str(data_archives[0]) + "/dyke.refactor.rk4.data")
RANGE_R = 0

try:
    RANGE_R = shel['RANGE_R']
finally:
    shel.close()

number_alive_global_start_np = np.zeros((RANGE_R, RANGE_R))
number_alive_local_start_np = np.zeros((RANGE_R, RANGE_R))
number_alive_global_end_np = np.zeros((RANGE_R, RANGE_R))
number_alive_local_end_np = np.zeros((RANGE_R, RANGE_R))
number_alive_start_np = np.zeros((RANGE_R, RANGE_R))
number_alive_end_np = np.zeros((RANGE_R, RANGE_R))

number_alive_global_start_ = []
number_alive_local_start_ = []
number_alive_global_end_ = []
number_alive_local_end_ = []
number_alive_start_ = []
number_alive_end_ = []
ENV_START_ = []


for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor.rk4.data")

    try:

        #RUN_ID = s['RUN_ID']

        number_alive_global_start = s['number_alive_global_start']
        number_alive_local_start = s['number_alive_local_start']
        number_alive_global_end = s['number_alive_global_end']
        number_alive_local_end = s['number_alive_local_end']
        number_alive_start = s['number_alive_start']
        number_alive_end = s['number_alive_end']

        ENV_START = s['ENV_START']
        RANGE_R = s['RANGE_R']

        number_alive_global_start_.append(number_alive_global_start)
        number_alive_local_start_.append(number_alive_local_start)
        number_alive_global_end_.append(number_alive_global_end)
        number_alive_local_end_.append(number_alive_local_end)
        number_alive_start_.append(number_alive_start)
        number_alive_end_.append(number_alive_end)
        ENV_START_.append(ENV_START)

        number_alive_global_start_np[ENV_START[0]][ENV_START[1]] += number_alive_global_start
        number_alive_local_start_np[ENV_START[0]][ENV_START[1]] += number_alive_local_start
        number_alive_global_end_np[ENV_START[0]][ENV_START[1]] += number_alive_global_end
        number_alive_local_end_np[ENV_START[0]][ENV_START[1]] += number_alive_local_end
        number_alive_start_np[ENV_START[0]][ENV_START[1]] += number_alive_start
        number_alive_end_np[ENV_START[0]][ENV_START[1]] += number_alive_end



        # SUPER DATA STRUCTURE NEEDED

        #(Eg, El, number_alive_global_start, number_alive_local_start, number_alive_start, number_alive_global_end, number_alive_local_end, number_alive_end)

        #Eg = [] ...
        #El = [] ...
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

def individual_plots():

    #===================================================================================================================
    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Global Alives Start')
    plt.xlabel('EL')
    plt.ylabel('EG')
    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
    idx = 0
    for start_vars in ENV_START_:
        heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_global_start_[idx]
        idx += 1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
    plt.show()
    #===================================================================================================================

    #===================================================================================================================
    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Local Alives Start')
    plt.xlabel('EL')
    plt.ylabel('EG')
    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
    idx = 0
    for start_vars in ENV_START_:
        idx += 1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
    plt.show()
    #===================================================================================================================

    #===================================================================================================================
    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Global + Local Alives Start')
    plt.xlabel('EL')
    plt.ylabel('EG')
    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
    idx = 0
    for start_vars in ENV_START_:
        heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_start_[idx]
        idx += 1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
    plt.show()
    #===================================================================================================================


    #===================================================================================================================
    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Global Alives End')
    plt.xlabel('EL')
    plt.ylabel('EG')
    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
    idx = 0
    for start_vars in ENV_START_:
        heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_global_end_[idx]
        idx += 1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
    plt.show()
    #===================================================================================================================

    #===================================================================================================================
    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Local Alives End')
    plt.xlabel('EL')
    plt.ylabel('EG')
    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
    idx = 0
    for start_vars in ENV_START_:
        heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_local_end_[idx]
        idx += 1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
    plt.show()
    #===================================================================================================================

    #===================================================================================================================
    plt.figure(figsize=(8,8), dpi=200)
    plt.title('Global + Local Alives End')
    plt.xlabel('EL')
    plt.ylabel('EG')
    heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
    idx = 0
    for start_vars in ENV_START_:
        heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_end_[idx]
        idx += 1
    plt.colorbar(plt.pcolor(heatmap))
    plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
    plt.show()
    #===================================================================================================================


if __name__ == '__main__':


    fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(20,10), dpi=300)
    plt.title('Global + Local Species Counts')

    minmin = np.min([np.min(number_alive_global_start_np),
                     np.min(number_alive_local_start_np),
                     np.min(number_alive_start_np),
                     np.min(number_alive_global_end_np),
                     np.min(number_alive_local_end_np),
                     np.min(number_alive_end_np)])
    maxmax = np.max([np.max(number_alive_global_start_np),
                     np.max(number_alive_local_start_np),
                     np.max(number_alive_start_np),
                     np.max(number_alive_global_end_np),
                     np.max(number_alive_local_end_np),
                     np.max(number_alive_end_np)])

    axes[0][0].title.set_text('Global Alive Start')
    axes[0][0].set_xlabel('El')
    axes[0][0].set_ylabel('Eg')
    axes[0][1].title.set_text('Local Alive Start')
    axes[0][1].set_xlabel('El')
    axes[0][1].set_ylabel('Eg')
    axes[0][2].title.set_text('Global + Local Alive Start')
    axes[0][2].set_xlabel('El')
    axes[0][2].set_ylabel('Eg')
    axes[1][0].title.set_text('Global Alive End')
    axes[1][0].set_xlabel('El')
    axes[1][0].set_ylabel('Eg')
    axes[1][1].title.set_text('Local Alive End')
    axes[1][1].set_xlabel('El')
    axes[1][1].set_ylabel('Eg')
    axes[1][2].title.set_text('Global + Local Alive End')
    axes[1][2].set_xlabel('El')
    axes[1][2].set_ylabel('Eg')

    im1 = axes[0][0].imshow(number_alive_global_start_np, vmin=minmin, vmax=maxmax,
                         extent=(0,100,0,100), aspect='equal', cmap='viridis')
    im2 = axes[0][1].imshow(number_alive_local_start_np, vmin=minmin, vmax=maxmax,
                         extent=(0,100,0,100), aspect='equal', cmap='viridis')
    im3 = axes[0][2].imshow(number_alive_start_np, vmin=minmin, vmax=maxmax,
                         extent=(0,100,0,100), aspect='equal', cmap='viridis')
    im4 = axes[1][0].imshow(number_alive_global_end_np, vmin=minmin, vmax=maxmax,
                         extent=(0,100,0,100), aspect='equal', cmap='viridis')
    im5 = axes[1][1].imshow(number_alive_local_end_np, vmin=minmin, vmax=maxmax,
                         extent=(0,100,0,100), aspect='equal', cmap='viridis')
    im6 = axes[1][2].imshow(number_alive_end_np, vmin=minmin, vmax=maxmax,
                         extent=(0,100,0,100), aspect='equal', cmap='viridis')

    # add space for colour bar
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    fig.colorbar(im6, cax=cbar_ax)
    plt.show()


    #fig, axes = plt.subplots(nrows=1, ncols=2)
    # generate randomly populated arrays
    #data1 = np.random.rand(10,10)*10
    #data2 = np.random.rand(10,10)*10 -7.5
    # find minimum of minima & maximum of maxima
    #minmin = np.min([np.min(data1), np.min(data2)])
    #maxmax = np.max([np.max(data1), np.max(data2)])
    #im1 = axes[0].imshow(data1, vmin=minmin, vmax=maxmax,
    #                     extent=(-5,5,-5,5), aspect='auto', cmap='viridis')
    #im2 = axes[1].imshow(data2, vmin=minmin, vmax=maxmax,
    #                     extent=(-5,5,-5,5), aspect='auto', cmap='viridis')
    # add space for colour bar
    #fig.subplots_adjust(right=0.85)
    #cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    #fig.colorbar(im2, cax=cbar_ax)
    #plt.show()

    #individual_plots()

    print("Completed")
