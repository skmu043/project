import sys, os
import shelve
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from tqdm import tqdm
import statistics
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase
import pandas as pd
from scipy import optimize
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


RESULT_DATA = []

data_dr = os.getcwd() + '/data_ztest_done'
data_archives = os.listdir(data_dr)

count = 0

plt.rcParams["font.family"] = "Times New Roman"

for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/ztest.exp.data")
    try:
        K_RESULT = s['K_RESULT']

        # K_RESULT.append((
        # 0 omega,
        # 1 mu,
        # 2 five,
        # 3 five_abundances,
        # 4 st,
        # 5 st_abundances,
        # 6 nw,
        # 7 nw_abundances
        # ))

        RESULT_DATA.append((K_RESULT[0][0],K_RESULT[0][1],K_RESULT[0][2],K_RESULT[0][3],K_RESULT[0][4],K_RESULT[0][5],K_RESULT[0][6],K_RESULT[0][7]))
        count += 1

    finally:
        s.close()


XFONT = 40
YFONT = 40
X_TICKS = 20
Y_TICKS = 20
XFIG = 30
YFIG = 30
TFONT = 40

def attractor_abundance():

    # K_RESULT.append((
    # 0 omega,
    # 1 mu,
    # 2 five,
    # 3 five_abundances,
    # 4 st,
    # 5 st_abundances,
    # 6 nw,
    # 7 nw_abundances
    # ))

    GAF_Attractor = []
    GAF_Abundance = []
    TGAF_Attractor = []
    TGAF_Abundance = []
    IGAF_Attractor = []
    IGAF_Abundance = []


    for data_point in RESULT_DATA:
        GAF  = data_point[3]
        for attractor_abundance in GAF:
            GAF_Attractor.append(attractor_abundance[0])
            GAF_Abundance.append(attractor_abundance[1])
        TGAF = data_point[5]
        for attractor_abundance in TGAF:
            TGAF_Attractor.append(attractor_abundance[0])
            TGAF_Abundance.append(attractor_abundance[1])
        IGAF = data_point[7]
        for attractor_abundance in IGAF:
            IGAF_Attractor.append(attractor_abundance[0])
            IGAF_Abundance.append(attractor_abundance[1])

    df1 = pd.DataFrame.from_dict(
        data=dict(GAF_Attractor=GAF_Attractor, GAF_Abundance=GAF_Abundance),
        orient='index',
    ).T
    df2 = pd.DataFrame.from_dict(
        data=dict(TGAF_Attractor=TGAF_Attractor, TGAF_Abundance=TGAF_Abundance),
        orient='index',
    ).T
    df3 = pd.DataFrame.from_dict(
        data=dict(IGAF_Attractor=IGAF_Attractor, IGAF_Abundance=IGAF_Abundance),
        orient='index',
    ).T

    plt.figure(figsize=(20,20), dpi=200)
    ax = plt.axes()
    plt.xlabel('Attractor', fontsize=40)
    plt.ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    #ax = sns.stripplot(data=df1, x="GAF_Attractor", y="GAF_Abundance")
    plt.scatter(GAF_Attractor, GAF_Abundance)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=6)
    plt.legend(prop={'size': 25},loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=5)
    plt.tight_layout()
    plt.savefig('GAF_attractor_abundance.jpg')
    plt.show()

    plt.figure(figsize=(20,20), dpi=200)
    ax = plt.axes()
    plt.xlabel('Attractor', fontsize=40)
    plt.ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    #ax = sns.stripplot(data=df2, x="TGAF_Attractor", y="TGAF_Abundance")
    plt.scatter(TGAF_Attractor, TGAF_Abundance)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=6)
    plt.legend(prop={'size': 25},loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=5)
    plt.tight_layout()
    plt.savefig('TGAF_attractor_abundance.jpg')
    plt.show()

    plt.figure(figsize=(20,20), dpi=200)
    ax = plt.axes()
    plt.xlabel('Attractor', fontsize=40)
    plt.ylabel('Total Abundance', fontsize=40)
    plt.xticks(fontsize=X_TICKS)
    plt.yticks(fontsize=Y_TICKS)
    #ax = sns.stripplot(data=df3, x="IGAF_Attractor", y="IGAF_Abundance")
    plt.scatter(IGAF_Attractor, IGAF_Abundance)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=8)
    ax.tick_params(which='minor', length=6)
    plt.legend(prop={'size': 25},loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=5)
    plt.tight_layout()
    plt.savefig('IGAF_attractor_abundance.jpg')
    plt.show()

attractor_abundance()


