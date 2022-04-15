import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)
shel = shelve.open(data_dr + "/" + str(data_archives[0]) + "/dyke.truncation.rk4.data")
RANGE_R = 0
times_step = 0

try:
    RANGE_R = shel['RANGE_R']
    times_steps = shel['times_steps']

finally:
    shel.close()

number_alive_global_start_np = np.zeros((RANGE_R, RANGE_R))
number_alive_global_end_np = np.zeros((RANGE_R, RANGE_R))
number_alive_start_np = np.zeros((RANGE_R, RANGE_R))
number_alive_end_np = np.zeros((RANGE_R, RANGE_R))
number_alive_diff_end_np = np.zeros((RANGE_R, RANGE_R))

number_alive_global_start_ = []
number_alive_global_end_ = []
number_alive_start_ = []
number_alive_end_ = []
ENV_START_ = []

results_gl = []
time_steps = []

for file in data_archives:
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.truncation.rk4.data")

    try:

        number_alive_global_start = s['number_alive_global_start']
        number_alive_global_end = s['number_alive_global_end']
        number_alive_start = s['number_alive_start']
        number_alive_end = s['number_alive_end']

        ENV_START = s['ENV_START']
        RANGE_R = s['RANGE_R']

        results = s['results']

        number_alive_global_start_.append(number_alive_global_start)
        number_alive_global_end_.append(number_alive_global_end)
        number_alive_start_.append(number_alive_start)
        number_alive_end_.append(number_alive_end)
        ENV_START_.append(ENV_START)

        results_gl.append((results[-2],results[-1]))

        np_eg = (RANGE_R) - (ENV_START[0])
        np_el = ENV_START[1] - 1
        if(ENV_START[0] == 0):
            np_eg = RANGE_R - 1
        if(ENV_START[1] == 0):
            np_el = ENV_START[1]
        if(ENV_START[0] == RANGE_R):
            np_eg = 0
        if(ENV_START[1] == RANGE_R):
            np_el = ENV_START[1] - 1


        number_alive_global_start_np[np_eg][np_el] += number_alive_global_start
        number_alive_global_end_np[np_eg][np_el] += number_alive_global_end
        number_alive_start_np[np_eg][np_el] += number_alive_start
        number_alive_end_np[np_eg][np_el] += number_alive_end

    finally:
        s.close()

number_alive_diff_end_np = number_alive_end_np - number_alive_start_np

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


if __name__ == '__main__':


    fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(10,10), dpi=300)
    plt.title('Global + Local Species (Start End) Difference')

    axes.title.set_text('Difference')
    axes.set_xlabel('El')
    axes.set_ylabel('Eg')

    minmin = abs(np.min([np.min(number_alive_diff_end_np)]))
    maxmax = abs(np.max([np.max(number_alive_diff_end_np)]))
    scale = maxmax
    if(minmin>=maxmax):     
        scale = minmin

    im1 = axes.imshow(number_alive_diff_end_np, extent=(0,RANGE_R,0,RANGE_R), aspect='equal', vmin=-scale, vmax=scale, cmap='Spectral')

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    fig.colorbar(im1, cax=cbar_ax)
    plt.savefig("3d_el_eg_alives_diff_" + str(random.randint(100, 999)) + ".png")
    plt.show()


    #===================================================================================================================

    #Verified Correct

    fig = plt.figure(figsize=(10,10), dpi=300)
    plt.title(label = "EL/EG")
    plt.xlabel('X - EL', fontsize=10)
    plt.ylabel('Y - EG', fontsize=10)

    for gl in results_gl:
        c_l = gl[1][-1]
        c_g = gl[0][-1]
        if(c_l < 0 or c_g < 0 or c_l > 100 or c_g > 100 ):
            plt.scatter(gl[1][0], gl[0][0],  marker='.', s=20, color=(float(0), float(0), float(1)))
            plt.plot(gl[1], gl[0],color=(float(0), float(0), float(1)))
        else:
            plt.scatter(gl[1][0], gl[0][0], marker='.', s=20, color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
            plt.plot(gl[1], gl[0], color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
            plt.scatter(gl[1][-1], gl[0][-1], marker='*', s=20, color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
    plt.savefig("2d_el_eg_" + str(random.randint(100, 999)) + ".png")
    plt.show()


    fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(20,10), dpi=300)
    plt.title('Global + Local Species Counts')

    minmin = np.min([np.min(number_alive_global_start_np),
                     np.min(number_alive_start_np),
                     np.min(number_alive_global_end_np),
                     np.min(number_alive_end_np)])
    maxmax = np.max([np.max(number_alive_global_start_np),
                     np.max(number_alive_start_np),
                     np.max(number_alive_global_end_np),
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
                         extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
    im2 = axes[0][1].imshow(number_alive_local_start_np, vmin=minmin, vmax=maxmax,
                         extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
    im3 = axes[0][2].imshow(number_alive_start_np, vmin=minmin, vmax=maxmax,
                         extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
    im4 = axes[1][0].imshow(number_alive_global_end_np, vmin=minmin, vmax=maxmax,
                         extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
    im5 = axes[1][1].imshow(number_alive_local_end_np, vmin=minmin, vmax=maxmax,
                         extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
    im6 = axes[1][2].imshow(number_alive_end_np, vmin=minmin, vmax=maxmax,
                         extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')

    # add space for colour bar
    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
    fig.colorbar(im6, cax=cbar_ax)
    plt.savefig("3d_el_eg_alives_" + str(random.randint(100, 999)) + ".png")
    plt.show()



    #===================================================================================================================


    print("Completed")
