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
from scipy.optimize import curve_fit
from _dbm import *
from scipy.ndimage.filters import gaussian_filter1d
from scipy.signal import savgol_filter


RESULT_DATA = []
UNIQ_SAMPLES = []

data_dr = os.getcwd() + '/data_attractors'
data_archives = os.listdir(data_dr)

#(K_list, omega,mu,five,seven,ten,st,nw)

for file in tqdm(data_archives):
    s = shelve.open(data_dr + "/" + str(file) + "/dyke.attractors.data")

    try:
        SAMPLE = s['K_RESULT']
        RESULT_DATA.append(SAMPLE)
    finally:
        s.close()

SPACE = 500

JI05 = [[] for _ in range(SPACE)]
JI07 = [[] for _ in range(SPACE)]
JI10 = [[] for _ in range(SPACE)]
ST05 = [[] for _ in range(SPACE)]
NW10 = [[] for _ in range(SPACE)]

for item in RESULT_DATA:
    for entry in item:
        #JI05[entry[0][0]-1].append(len(entry[3]))
        #JI07[entry[0][0]-1].append(len(entry[4]))
        #JI10[entry[0][0]-1].append(len(entry[5]))
        ST05[entry[0][0]-1].append(len(entry[6]))
        NW10[entry[0][0]-1].append(len(entry[7]))

biotic_components = []
number_of_attractors = []
specific_model = []

for _ in range(SPACE):
    #for item in JI05[_]:
    #    biotic_components.append(_+1)
    #    number_of_attractors.append(item)
    #    specific_model.append("JI05")
        #for item in JI07[_]:
    #    biotic_components.append(_+1)
    #    number_of_attractors.append(item)
    #    specific_model.append("JI07")
    #for item in JI10[_]:
    #    biotic_components.append(_+1)
    #    number_of_attractors.append(item)
    #    specific_model.append("JI10")

    #for item in ST05[_]:
    #    biotic_components.append(_+1)
    #    number_of_attractors.append(item)
    #    specific_model.append("ST05")
    for item in NW10[_]:
        biotic_components.append(_+1)
        number_of_attractors.append(item)
        specific_model.append("NW10")



zipped = list(zip(biotic_components,number_of_attractors,specific_model))
df = pd.DataFrame(zipped, columns=['biotic_components','number_of_attractors','specific_model'])

#print(df.groupby(['biotic_components']).mean())
#print(df.groupby(['specific_model', 'biotic_components']).std())
df2 = df.groupby(['biotic_components']).mean()
df3 = df.groupby(['biotic_components']).std()

print(df3.keys())

df3_t = df2['number_of_attractors'] - df3['number_of_attractors']
df3_b = df2['number_of_attractors'] + df3['number_of_attractors']

#print(df2.keys())

#print(df2.index)
fig, ax = plt.subplots(figsize=(20,20), dpi= 200)

x = df2.index.tolist()
y = df2['number_of_attractors'].tolist()
#ax.set_xscale('log')
#plt.plot(x,y)
yhat = savgol_filter(y, 51, 3)

#plt.plot(x,yhat_)
sns.lineplot(x=x,y=y)
sns.lineplot(x=x, y=yhat)

fill1 = savgol_filter(df3_t.tolist(), 51, 3)
fill2 = savgol_filter(df3_b.tolist(), 51, 3)
ax.fill_between(x, fill1, fill2, color="blue", alpha=0.3)
sns.lineplot(x=x, y=fill1)
sns.lineplot(x=x, y=fill2)

ax.set_xscale('log')
plt.show()

#Final_array_smooth = gaussian_filter1d(df2["number_of_attractors"], sigma=2)

#fig, ax = plt.subplots(figsize=(20,20), dpi= 200)
#sns.stripplot(x="biotic_components", y = 'number_of_attractors', hue='specific_model', data=df, jitter=0.25)
#sns.pointplot(data = df, x = 'biotic_components', y = 'number_of_attractors', hue='specific_model',  dodge=True, join=True, errwidth = 3, capsize = 0.5, markersize = 50)
#sns.lmplot(data = df, x = 'biotic_components', y = 'number_of_attractors', hue='specific_model');
#sns.lineplot(data = df, x = 'biotic_components', y = 'number_of_attractors', hue='specific_model');
#sns.regplot(data = df, x = 'biotic_components', y = 'number_of_attractors', scatter = False)

#ax.set_xscale('log')
#fig.update_layout(xaxis_type="log")
#ax.set_xticklabels(['1','2', '5', '10'])
#plt.tight_layout()
#plt.savefig('attractors_JI05_JI07_JI10_ST05_NW10.jpg' )
#plt.show()



