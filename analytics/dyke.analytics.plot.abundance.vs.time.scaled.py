import shelve, os
import matplotlib.pyplot as plt

#plot abundance of each species over time where abundance is scaled up by R
def plot_aot_scaled():
    plt.figure(figsize=(30,30))
    plt.title('Abundance + Temp over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance Scaled UP + Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-50, R+20)
    plt.xlim(0, end)
    for x in range(K):
        plt.plot(time,rAxR[x],label = 'id %s'%x)

    plt.plot(time,rE, 'r.', label='E')
    plt.show()

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):

    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke.data.db")

    try :
        rAxR            = s['rAxR']
        time            = s['time']
        rE              = s['rE']
        K               = s['K']
        R               = s['R']
        end             = s['end']

        plot_aot_scaled()      #plot abundance of each species over time scaled by R

    finally:
        s.close()

else:
    print("Invalid Selection")



