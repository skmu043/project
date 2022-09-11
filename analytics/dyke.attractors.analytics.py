import os
import shelve
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
import pandas as pd
from scipy.signal import savgol_filter
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

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

biotic_components = []
number_of_attractors = []
specific_model = []

for item in RESULT_DATA:
    for entry in item:

        biotic_components.append(entry[0][0]-1)
        number_of_attractors.append(len(entry[3]))
        specific_model.append("Niche = 5")

        biotic_components.append(entry[0][0]-1)
        number_of_attractors.append(len(entry[4]))
        specific_model.append("Niche = 7")

        biotic_components.append(entry[0][0]-1)
        number_of_attractors.append(len(entry[5]))
        specific_model.append("Niche = 10")

        #biotic_components.append(entry[0][0]-1)
        #number_of_attractors.append(len(entry[6]))
        #specific_model.append("ST05")

        #biotic_components.append(entry[0][0]-1)
        #number_of_attractors.append(len(entry[7]))
        #specific_model.append("NW10")


zipped = list(zip(biotic_components,number_of_attractors,specific_model))
df = pd.DataFrame(zipped, columns=['biotic_components','number_of_attractors','specific_model'])

#JI05 = df[df['specific_model'] == 'JI05']
#JI07 = df[df['specific_model'] == 'JI07']
#JI10 = df[df['specific_model'] == 'JI10']
#ST05 = df[df['specific_model'] == 'ST05']
#NW10 = df[df['specific_model'] == 'NW10']

fig, ax = plt.subplots(figsize=(20,20), dpi= 200)


for specific_model in df['specific_model'].unique():

    SMDF_mean = df[df['specific_model'] == specific_model].groupby(['biotic_components']).mean()
    SMDF_std = df[df['specific_model'] == specific_model].groupby(['biotic_components']).std()

    SMDF_upper = SMDF_mean['number_of_attractors'] - SMDF_std['number_of_attractors']
    SMDF_lower = SMDF_mean['number_of_attractors'] + SMDF_std['number_of_attractors']

    x = SMDF_mean.index.tolist()
    y = SMDF_mean['number_of_attractors'].tolist()

    yhat = savgol_filter(y, 51, 3)
    yhat[0]=0
    fill1 = savgol_filter(SMDF_upper.tolist(), 51, 3)
    fill2 = savgol_filter(SMDF_lower.tolist(), 51, 3)

    #sns.lineplot(x=x,y=y, linewidth=3)
    #sns.scatterplot(x=x,y=y)
    #plt.plot(x,y)
    #ax.fill_between(x, fill1, fill2, alpha=0.3)
    sns.lineplot(x=x, y=y, label = str(specific_model), markers=["o", "o", "o"])
    #sns.lineplot(x=x, y=fill1)
    #sns.lineplot(x=x, y=fill2)

#ax.set_xscale('log')

plt.title('Attractors for three niche sizes', fontsize=40)
ax.set_xlabel('Number of Species', fontsize=40)
ax.set_ylabel('Number of Attractors', fontsize=40)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.tick_params(which='both', width=1)
ax.tick_params(which='major', length=7)
ax.tick_params(which='minor', length=4)
ax.set_xscale('log')

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(loc=2,prop={'size': 20})
plt.tight_layout()
plt.savefig('attractors_JI05_JI07_JI10_ST05_NW10.jpg')
plt.show()


#df2 = df.groupby(['biotic_components']).mean()
#df3 = df.groupby(['biotic_components']).std()

#df3_t = df2['number_of_attractors'] - df3['number_of_attractors']
#df3_b = df2['number_of_attractors'] + df3['number_of_attractors']

#fig, ax = plt.subplots(figsize=(20,20), dpi= 200)

#x = df2.index.tolist()
#y = df2['number_of_attractors'].tolist()

#yhat = savgol_filter(y, 51, 3)
#yhat = savgol_filter(yhat, 51, 3)
#fill1 = savgol_filter(df3_t.tolist(), 51, 3)
#fill1 = savgol_filter(fill1, 51, 3)
#fill2 = savgol_filter(df3_b.tolist(), 51, 3)
#fill2 = savgol_filter(fill2, 51, 3)
#sns.lineplot(x=x,y=y)
#ax.fill_between(x, fill1, fill2, alpha=0.3)
#sns.lineplot(x=x, y=yhat)
#sns.lineplot(x=x, y=fill1
#sns.lineplot(x=x, y=fill2)

#ax.set_xscale('log')
#plt.show()
