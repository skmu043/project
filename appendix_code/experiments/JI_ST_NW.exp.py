import sys
import shelve
import math
import numpy as np

# This file gets called with input parameters for each experiment
# Checks input parameters
if int(len(sys.argv)) != int(2):
    print("Args: shelve file name which contains all of > (K, R, P, E, start, end, step, EN, OE, RUN_ID)")
    print("e.g K=100, R=100, P=0, E=X, start=0, end=200, step=0.01, EN=2, OE=5, RUN_ID : epoch")
    print("exit")
    sys.exit()

s = shelve.open(str(sys.argv[1]))

try:
    # Read experiment starting parameters
    SPECIES_K           = s['SPECIES_K']
    RANGE_R             = s['RANGE_R']
    TIME_START          = s['TIME_START']
    TIME_END            = s['TIME_END']
    TIME_STEP           = s['TIME_STEP']
    ENV_VARS            = s['ENV_VARS']
    NICHE               = s['NICHE']
    SURVIVAL_THRESHOLD  = s['SURVIVAL_THRESHOLD']
    ENV_START           = s['ENV_START']
    omega               = s['omega']
    mu                  = s['mu']
    exp_name            = s['exp_name']
    data_directory      = s['data_directory']
    shelve_file         = s['shelve_file']

finally:
    s.close()

# Initilization of the system state that is used for the RK4 simulation
system_state        = np.zeros(SPECIES_K+ENV_VARS)
SystemTemperature   = ENV_START

# ---------------------------------------------------------------------------------------------------------------------#
# NICHE - Match Range on Niche Change

def fYaI(Xe, Ni, u, T):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if(abundance <= T):
        abundance = 0

    return(abundance)

def temperature_range(Ya, Ni, u):
    return (u + (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))))

def negative_temperature_range(Ya, Ni, u):
    return (u - (math.sqrt(((math.log(Ya,math.e) / -1) * (2*(Ni**2))))))

def abundance_temperature_range(Xe, Ni, u, NRange):

    abundance = ((math.e) ** ((-1) * (((abs(Xe-u)) ** 2) / (2*(Ni**2)))))

    if((Xe >= u + NRange) or (Xe <= u - NRange)):
        abundance = 0

    return(abundance)
# ---------------------------------------------------------------------------------------------------------------------#
# Initialization of the abundance values needed for the first iteration
# Processes abundance calculation for the JI and ST model
if(NICHE == 5):
    for s_i in range(SPECIES_K):
        a_star = np.exp(- abs(SystemTemperature-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
        if a_star < SURVIVAL_THRESHOLD:
            a_star = 0
        system_state[s_i] = a_star

# Processes abundance calculation for the NW
if(NICHE == 10):
    for s_i in range(SPECIES_K):
        # Returns the distance from the optimal growing temperature to the edge of the abundance graph
        NRange = (temperature_range(SURVIVAL_THRESHOLD, 5, mu[0][s_i]) - mu[0][s_i])
        # Utilizes the distance to calculate the abundance for a species
        a_star = abundance_temperature_range(SystemTemperature, NICHE, mu[0][s_i], NRange)
        system_state[s_i] = a_star
# ---------------------------------------------------------------------------------------------------------------------#
# Captures the total starting abundance that is calculated above
abundance_start = sum(system_state)

# Determines the number of alive species at the start of the simulation
number_alive_start = 0
for abundance_init in system_state:
    if abundance_init > 0:
        number_alive_start +=1

# Initializes system temperature in the system state
for _ in range(ENV_VARS):
    system_state[SPECIES_K+_] = ENV_START

def rates_of_change_system_state(system_state):

    # This function updates all the system state variables
    rate_of_change = system_state.copy()

    SystemTemperature = system_state[SPECIES_K+0]
    # Abundance Update for the JI and ST model
    if(NICHE == 5):
        for s_i in range(SPECIES_K):
            a_star = np.exp(- abs(SystemTemperature-mu[0][s_i]) ** 2 / ( 2 * NICHE ** 2 ))
            if a_star < SURVIVAL_THRESHOLD:
                a_star = 0
            rate_of_change[s_i] =  a_star - system_state[s_i]
    # Abundance Update for NW model
    if(NICHE == 10):
        for s_i in range(SPECIES_K):
            NRange = (temperature_range(SURVIVAL_THRESHOLD, 5, mu[0][s_i]) - mu[0][s_i])
            a_star = abundance_temperature_range(SystemTemperature, NICHE, mu[0][s_i], NRange)
            rate_of_change[s_i] =  a_star - system_state[s_i]

    biotic_force_FG = 0

    # Calculates the biotic force the species are applying to the system temperature
    for s_i in range(SPECIES_K):
        biotic_force_FG += (system_state[s_i] * omega[0][s_i])

    rate_of_change[SPECIES_K+0] = (biotic_force_FG)

    return(rate_of_change)

if __name__ == '__main__':

    results = [[] for _ in range(SPECIES_K+ENV_VARS)]
    times_steps=[]

    for step in np.arange(TIME_START, TIME_END, TIME_STEP):
        times_steps.append(step)
        for _ in range(SPECIES_K+ENV_VARS):
            results[_].append(system_state[_])

        # Uses the Runge–Kutta method to simulate the system for 200 time steps
        k1 = TIME_STEP * rates_of_change_system_state(system_state)
        k2 = TIME_STEP * rates_of_change_system_state(system_state + k1 * 0.5)
        k3 = TIME_STEP * rates_of_change_system_state(system_state + k2 * 0.5)
        k4 = TIME_STEP * rates_of_change_system_state(system_state + k3)

        system_state += ((k1 + (2*k2) + (2*k3) + k4)/6)

    number_alive_end = 0
    aliveness = 0
    abundance_end = 0
    # Updates the number alive and total abundance at the end of the simulation
    for abundance_stream in results[:-1]:
        if abundance_stream[-1] > aliveness:
            number_alive_end +=1
            abundance_end += abundance_stream[-1]


    s = shelve.open(str(sys.argv[1]))

    try:
        # Write the experiment results to the shelve file for analysis
        s['NUMBER_ALIVE_START']     = number_alive_start            # --- Number of alive species at the start
        s['NUMBER_ALIVE_END']       = number_alive_end              # --- Number of alive species at the end
        s['TOTAL_ABUNDANCE_START']  = abundance_start               # --- Total abundance at the start
        s['TOTAL_ABUNDANCE_END']    = abundance_end                 # --- Total abundance at the end
        s['ENV_END']                = results[SPECIES_K][-1]        # --- System Temperature at the end

    finally:
        s.close()