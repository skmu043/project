import shelve, os
import matplotlib.pyplot as plt

#plot species that increase temperature and decrease temperature
def plot_aot_inc_dec():
    plt.figure(figsize=(20,10))
    plt.title('Species Abundance (Blues Decrease Temperature while Reds Increase)', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for x in range(K):
        if(w[x]==0):
            plt.plot(time,rAx[x],'k-')
        if(w[x]<0):
            plt.plot(time,rAx[x],'b-')
        if(w[x]>0):
            plt.plot(time,rAx[x],'r-')
    plt.show()

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):
    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke.data.db")

    try :

        args            = s['sys.argv']
        temperatures    = s['temperatures']
        biotic_force    = s['biotic_force']
        w               = s['w']
        u               = s['u']
        rAxR            = s['rAxR']
        time            = s['time']
        rE              = s['rE']
        rF              = s['rF']
        rP              = s['rP']
        rEt             = s['rEt']
        rAx             = s['rAx']

        K               = s['K']
        R               = s['R']
        P               = s['P']
        E               = s['E']
        start           = s['start']
        end             = s['end']
        step            = s['step']

        plot_aot_inc_dec()      #plot species that increase temperature and decrease temperature

    finally:
        s.close()

else:
    print("Invalid Selection")



