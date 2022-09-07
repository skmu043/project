import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm

data_dr = os.getcwd() + '/data_global_local'
data_archives = os.listdir(data_dr)
#shel = shelve.open(data_dr + "/" + str(data_archives[0]) + "/dyke.refactor.rk4.data")

SPECIES_K = 0
RANGE_R = 0
omega = []
mu = []
local_population_index = []
ENS_START = []
results = []

DATA_RESULTS = [] # omega, mu, (JI_ST), local_population_index, results

for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.refactor.rk4.data")

    try:

        SPECIES_K = s['SPECIES_K']
        RANGE_R = s['RANGE_R']
        omega = s['omega']
        mu = s['mu']
        local_population_index = s['local_population_index']
        ENS_START = s['ENV_START']
        JI_ST = s['JI_ST']
        results = s['results']

        DATA_RESULTS.append((omega, mu, JI_ST, local_population_index, results))

    finally:
        s.close()
#
# [0] omega
# [1] mu
# [2] ji_st
# [3] local_population_index
# [4] abundance of each of the k species with E1 and E2
#  --- [0] ... [99]
#         [K1 abundance]
#  --- [100] .. [101]
#          [E1]
#          [E2]


UNIQUE_SAMPLES = []

for data_entry in DATA_RESULTS:
    if((data_entry[0], data_entry[1]) not in UNIQUE_SAMPLES):
        UNIQUE_SAMPLES.append((data_entry[0], data_entry[1]))

print(len(UNIQUE_SAMPLES))

for data_entry in DATA_RESULTS:


    x = np.linspace(0, 200, 200)
    #print(x)

    r1 = data_entry[4]

    index = 0
    for item in r1:
        #print(item)
        if(index <= 99):
            plt.plot(x, item)

        index+=1
    plt.savefig(str(index) + ".jpg")
    #plt.show()

#
# def individual_plots():
#
#     #===================================================================================================================
#     plt.figure(figsize=(8,8), dpi=200)
#     plt.title('Global Alives Start')
#     plt.xlabel('EL')
#     plt.ylabel('EG')
#     heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
#     idx = 0
#     for start_vars in tqdm(ENV_START_):
#         heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_global_start_[idx]
#         idx += 1
#     plt.colorbar(plt.pcolor(heatmap))
#     plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
#     plt.show()
#     #===================================================================================================================
#
#     #===================================================================================================================
#     plt.figure(figsize=(8,8), dpi=200)
#     plt.title('Local Alives Start')
#     plt.xlabel('EL')
#     plt.ylabel('EG')
#     heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
#     idx = 0
#     for start_vars in tqdm(ENV_START_):
#         idx += 1
#     plt.colorbar(plt.pcolor(heatmap))
#     plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
#     plt.show()
#     #===================================================================================================================
#
#     #===================================================================================================================
#     plt.figure(figsize=(8,8), dpi=200)
#     plt.title('Global + Local Alives Start')
#     plt.xlabel('EL')
#     plt.ylabel('EG')
#     heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
#     idx = 0
#     for start_vars in tqdm(ENV_START_):
#         heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_start_[idx]
#         idx += 1
#     plt.colorbar(plt.pcolor(heatmap))
#     plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
#     plt.show()
#     #===================================================================================================================
#
#
#     #===================================================================================================================
#     plt.figure(figsize=(8,8), dpi=200)
#     plt.title('Global Alives End')
#     plt.xlabel('EL')
#     plt.ylabel('EG')
#     heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
#     idx = 0
#     for start_vars in tqdm(ENV_START_):
#         heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_global_end_[idx]
#         idx += 1
#     plt.colorbar(plt.pcolor(heatmap))
#     plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
#     plt.show()
#     #===================================================================================================================
#
#     #===================================================================================================================
#     plt.figure(figsize=(8,8), dpi=200)
#     plt.title('Local Alives End')
#     plt.xlabel('EL')
#     plt.ylabel('EG')
#     heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
#     idx = 0
#     for start_vars in tqdm(ENV_START_):
#         heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_local_end_[idx]
#         idx += 1
#     plt.colorbar(plt.pcolor(heatmap))
#     plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
#     plt.show()
#     #===================================================================================================================
#
#     #===================================================================================================================
#     plt.figure(figsize=(8,8), dpi=200)
#     plt.title('Global + Local Alives End')
#     plt.xlabel('EL')
#     plt.ylabel('EG')
#     heatmap = [[0 for _ in range(RANGE_R)] for _ in range(RANGE_R)]
#     idx = 0
#     for start_vars in tqdm(ENV_START_):
#         heatmap[int(start_vars[0])][int(start_vars[1])] += number_alive_end_[idx]
#         idx += 1
#     plt.colorbar(plt.pcolor(heatmap))
#     plt.imshow(heatmap, cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=12)
#     plt.show()
#     #===================================================================================================================
#
#
# if __name__ == '__main__':
#
#
#     fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(10,10), dpi=300)
#     plt.title('Global + Local Species (Start End) Difference')
#
#     axes.title.set_text('Difference')
#     axes.set_xlabel('El')
#     axes.set_ylabel('Eg')
#
#     minmin = abs(np.min([np.min(number_alive_diff_end_np)]))
#     maxmax = abs(np.max([np.max(number_alive_diff_end_np)]))
#     scale = maxmax
#     if(minmin>=maxmax):
#         scale = minmin
#
#     im1 = axes.imshow(number_alive_diff_end_np, extent=(0,RANGE_R,0,RANGE_R), aspect='equal', vmin=-scale, vmax=scale, cmap='Spectral')
#
#     fig.subplots_adjust(right=0.85)
#     cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
#     fig.colorbar(im1, cax=cbar_ax)
#     plt.savefig("3d_el_eg_alives_diff_" + str(random.randint(100, 999)) + ".png")
#     plt.show()
#
#
#     #===================================================================================================================
#
#     #Verified Correct
#
#     fig = plt.figure(figsize=(10,10), dpi=300)
#     plt.title(label = "EL/EG")
#     plt.xlabel('X - EL', fontsize=10)
#     plt.ylabel('Y - EG', fontsize=10)
#
#     for gl in tqdm(results_gl):
#         c_l = gl[1][-1]
#         c_g = gl[0][-1]
#         if(c_l < 0 or c_g < 0 or c_l > 100 or c_g > 100 ):
#             plt.scatter(gl[1][0], gl[0][0],  marker='.', s=20, color=(float(0), float(0), float(1)))
#             plt.plot(gl[1], gl[0],color=(float(0), float(0), float(1)))
#         else:
#             plt.scatter(gl[1][0], gl[0][0], marker='.', s=20, color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
#             plt.plot(gl[1], gl[0], color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
#             plt.scatter(gl[1][-1], gl[0][-1], marker='*', s=20, color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
#     plt.savefig("2d_el_eg_" + str(random.randint(100, 999)) + ".png")
#     plt.show()
#
#     #===================================================================================================================
#
#     #fig = plt.figure(figsize=(10,10), dpi=300)
#     #ax = fig.add_subplot(111, projection='3d', adjustable='box')
#     #ax.set_title(label = "EL/EG over Time Steps")
#     #ax.set_xlabel('X - EL', fontsize=10)
#     #ax.set_ylabel('Y - EG', fontsize=10)
#     #ax.set_zlabel('Z - Time Steps', fontsize=10)
#     #ax.set_xlim([-10,RANGE_R])
#     #ax.set_ylim([-10,RANGE_R])
#
#     #for gl in results_gl:
#     #    c_l = gl[1][-1]
#     #    c_g = gl[0][-1]
#     #    if(c_l < 0 or c_g < 0 or c_l > 100 or c_g > 100 ):
#     #        ax.scatter(gl[1][0], gl[0][0], times_steps[0], marker='.', s=20, color=(float(0), float(0), float(1)))
#     #        ax.plot(gl[1], gl[0], times_steps,color=(float(0), float(0), float(1)))
#     #    else:
#     #        ax.scatter(gl[1][0], gl[0][0], times_steps[0], marker='.', s=20, color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
#     #        ax.plot(gl[1], gl[0], times_steps,color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
#     #        ax.scatter(gl[1][-1], gl[0][-1], times_steps[-1], marker='*', s=20, color=(float(gl[1][-1]/100), float(gl[0][-1]/100), float(0.5)))
#     #plt.savefig("3d_el_eg_" + str(random.randint(100, 999)) + ".png")
#     #plt.show()
#
#     #===================================================================================================================
#
#
#     fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(20,10), dpi=300)
#     plt.title('Global + Local Species Counts')
#
#     minmin = np.min([np.min(number_alive_global_start_np),
#                      np.min(number_alive_local_start_np),
#                      np.min(number_alive_start_np),
#                      np.min(number_alive_global_end_np),
#                      np.min(number_alive_local_end_np),
#                      np.min(number_alive_end_np)])
#     maxmax = np.max([np.max(number_alive_global_start_np),
#                      np.max(number_alive_local_start_np),
#                      np.max(number_alive_start_np),
#                      np.max(number_alive_global_end_np),
#                      np.max(number_alive_local_end_np),
#                      np.max(number_alive_end_np)])
#
#     axes[0][0].title.set_text('Global Alive Start')
#     axes[0][0].set_xlabel('El')
#     axes[0][0].set_ylabel('Eg')
#     axes[0][1].title.set_text('Local Alive Start')
#     axes[0][1].set_xlabel('El')
#     axes[0][1].set_ylabel('Eg')
#     axes[0][2].title.set_text('Global + Local Alive Start')
#     axes[0][2].set_xlabel('El')
#     axes[0][2].set_ylabel('Eg')
#     axes[1][0].title.set_text('Global Alive End')
#     axes[1][0].set_xlabel('El')
#     axes[1][0].set_ylabel('Eg')
#     axes[1][1].title.set_text('Local Alive End')
#     axes[1][1].set_xlabel('El')
#     axes[1][1].set_ylabel('Eg')
#     axes[1][2].title.set_text('Global + Local Alive End')
#     axes[1][2].set_xlabel('El')
#     axes[1][2].set_ylabel('Eg')
#
#     im1 = axes[0][0].imshow(number_alive_global_start_np, vmin=minmin, vmax=maxmax,
#                          extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
#     im2 = axes[0][1].imshow(number_alive_local_start_np, vmin=minmin, vmax=maxmax,
#                          extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
#     im3 = axes[0][2].imshow(number_alive_start_np, vmin=minmin, vmax=maxmax,
#                          extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
#     im4 = axes[1][0].imshow(number_alive_global_end_np, vmin=minmin, vmax=maxmax,
#                          extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
#     im5 = axes[1][1].imshow(number_alive_local_end_np, vmin=minmin, vmax=maxmax,
#                          extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
#     im6 = axes[1][2].imshow(number_alive_end_np, vmin=minmin, vmax=maxmax,
#                          extent=(0,RANGE_R,0,RANGE_R), aspect='equal', cmap='viridis')
#
#     # add space for colour bar
#     fig.subplots_adjust(right=0.85)
#     cbar_ax = fig.add_axes([0.88, 0.15, 0.04, 0.7])
#     fig.colorbar(im6, cax=cbar_ax)
#     plt.savefig("3d_el_eg_alives_" + str(random.randint(100, 999)) + ".png")
#     plt.show()
#
#     individual_plots()
#
#
#     #===================================================================================================================
#
#
#     print("Completed")
