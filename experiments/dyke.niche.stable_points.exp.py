import sys
import shelve
import math
import numpy as np
from scipy import optimize

s = shelve.open(str(sys.argv[1]))
omega = []
mu = []

SPECIES_K   = 100                  # ----------- Number of Biotic Components
RANGE_R     = 100                  # ----------- Essential Range
TIME_START  = 0                     # ----------- Start of Simulation
TIME_END    = 200                   # ----------- Length of Simulation
TIME_STEP   = 1                   # ----------- Time Step3
ENV_VARS    = 1                     # ----------- Number of Environment Variables
NICHE = 5                           # ----------- Niche Size
LOCAL_SIZE  = 50                    # ----------- Local Population Size (%)
ALIVE_THRESHOLD = 0
ENV_START=[0]


try:
    exp_name            = s['exp_name']
    omega            = s['omega']
    mu            = s['mu']

finally:
    s.close()


########################################################################################################################
def f1(x):
    biotic_force = []
    for y in range(SPECIES_K):
        biotic_force.append(((math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))) * omega[0][y])
    return(np.sum((np.array(biotic_force, dtype=float))))

def ji_stable_point_return():


    temp_range = np.arange(-25, 125, 0.01)
    biotic_values = []
    for temp in temp_range:
        biotic_values.append(f1(temp))

    fsolve_search_points = []
    positive_or_negative_current = 0
    # get first point ...
    first_point = 0
    index_first_point = 0
    for biotic_point in biotic_values:
        if(biotic_point>0):
            first_point = biotic_point
            break
        if(biotic_point<0):
            first_point = biotic_point
            break
        index_first_point+=1

    if(first_point > 0):
        positive_or_negative_current = 1
    if(first_point < 0):
        positive_or_negative_current = -1


    for check_index in range(index_first_point+1, len(biotic_values)):

        if(biotic_values[check_index]>0):
            if(positive_or_negative_current == -1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = +1
        if(biotic_values[check_index]<0):
            if(positive_or_negative_current == 1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = -1

    fsolve_final_roots = []
    for crossing_zero in fsolve_search_points:
        roots = optimize.fsolve(f1,crossing_zero)
        fsolve_final_roots.append(roots)

    final_stable_points = []
    for each_point in fsolve_final_roots:
        if(f1(each_point-0.01) > 0 and f1(each_point+0.01) < 0):
            final_stable_points.append(each_point)

    biotic_force = [[] for _ in range(SPECIES_K)]
    step = 0.01

    results_points = []
    for stable_p in final_stable_points:
        results_points.append(stable_p[0])

    return(set(results_points))

def f1x(x):
    biotic_force = []
    truncation=0.2

    for y in range(SPECIES_K):
        aliveness = (math.e) ** ((-1) * (((abs(x-mu[0][y])) ** 2) / (2*(NICHE**2))))
        if(aliveness <= truncation and aliveness >= (-1 * truncation)):
            biotic_force.append(0)
        else:
            biotic_force.append(aliveness * omega[0][y])

    return(np.sum((np.array(biotic_force, dtype=float))))

def st_stable_point_return():

    temp_range = np.arange(-25, 125, 0.01)
    biotic_values = []
    for temp in temp_range:
        biotic_values.append(f1x(temp))

    fsolve_search_points = []
    positive_or_negative_current = 0
    # get first point ...
    first_point = 0
    index_first_point = 0
    for biotic_point in biotic_values:
        if(biotic_point>0):
            first_point = biotic_point
            break
        if(biotic_point<0):
            first_point = biotic_point
            break
        index_first_point+=1

    if(first_point > 0):
        positive_or_negative_current = 1
    if(first_point < 0):
        positive_or_negative_current = -1


    for check_index in range(index_first_point+1, len(biotic_values)):

        if(biotic_values[check_index]>0):
            if(positive_or_negative_current == -1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = +1
        if(biotic_values[check_index]<0):
            if(positive_or_negative_current == 1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = -1

    fsolve_final_roots = []
    for crossing_zero in fsolve_search_points:
        roots = optimize.fsolve(f1x,crossing_zero)
        fsolve_final_roots.append(roots)

    final_stable_points = []
    for each_point in fsolve_final_roots:
        if(f1x(each_point-0.01) > 0 and f1x(each_point+0.01) < 0):
            final_stable_points.append(each_point)


    results_points = []
    for stable_p in final_stable_points:
        results_points.append(stable_p[0])
    return(set(results_points))

########################################################################################################################

def fYa(Xe, Ni, u):
    return (((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2))))))

def fXe(Ya, Ni, u):
    return (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))) + u

def fYaI(Xe, Ni, u, T):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if(abundance <= T):
        abundance = 0

    return(abundance)

def fXe(Ya, Ni, u):
    return (u + (math.sqrt(((-1 * math.log(Ya,math.e)) * (2*(Ni**2))))))

def fXe_negative(Ya, Ni, u):
    return (u - (math.sqrt(((-1 * math.log(Ya,math.e)) * (2*(Ni**2))))))

def fYaIx(Xe, Ni, u, NRange):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if((Xe >= u + NRange) or (Xe <= u - NRange)):
        abundance = 0

    return(abundance)
########################################################################################################################

def f1x2(x):

    biotic_force = []
    truncation=0.2

    for y in range(SPECIES_K):
        NRange = (fXe(0.2, 5, mu[0][y]) - mu[0][y])
        aliveness = fYaIx(x, 10, mu[0][y], NRange)
        biotic_force.append(aliveness * omega[0][y])

    return(np.sum((np.array(biotic_force, dtype=float))))

def nw_stable_point_return():


    temp_range = np.arange(-25, 125, 0.01)
    biotic_values = []
    for temp in temp_range:
        biotic_values.append(f1x2(temp))

    fsolve_search_points = []
    positive_or_negative_current = 0
    # get first point ...
    first_point = 0
    index_first_point = 0
    for biotic_point in biotic_values:
        if(biotic_point>0):
            first_point = biotic_point
            break
        if(biotic_point<0):
            first_point = biotic_point
            break
        index_first_point+=1

    if(first_point > 0):
        positive_or_negative_current = 1
    if(first_point < 0):
        positive_or_negative_current = -1


    for check_index in range(index_first_point+1, len(biotic_values)):

        if(biotic_values[check_index]>0):
            if(positive_or_negative_current == -1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = +1
        if(biotic_values[check_index]<0):
            if(positive_or_negative_current == 1): # sign change has happened
                fsolve_search_points.append(temp_range[check_index])
                positive_or_negative_current = -1

    reduced_points = []
    for each_point in fsolve_search_points:
        if(f1x2(each_point-0.01) > 0 and f1x2(each_point+0.01) < 0):
            reduced_points.append(each_point)

    fsolve_final_roots = []
    for crossing_zero in reduced_points:
        roots = optimize.fsolve(f1x2,crossing_zero)
        fsolve_final_roots.append(roots)

    final_stable_points = reduced_points

    return(set(final_stable_points))




if __name__ == '__main__':

    ji = ji_stable_point_return()
    st = st_stable_point_return()
    nw = nw_stable_point_return()


    s = shelve.open(str(sys.argv[1]))


    try:
        s['JI_STABLE_POINTS'] = ji
        s['ST_STABLE_POINTS'] = st
        s['NW_STABLE_POINTS'] = nw

    finally:
        s.close()


