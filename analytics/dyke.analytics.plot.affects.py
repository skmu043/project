import shelve, os
import matplotlib.pyplot as plt

#plot affects values for each species
def plot_w():

    plt.figure(figsize=(20,5))
    plt.ylim(-1, 1)
    plt.title('Affects values for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Affects', fontsize=18)
    plt.plot(w, 'k.', label='w')
    plt.legend(loc=5, prop={'size': 30})
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

        plot_w()               #plot affects values for each species

    finally:
        s.close()

else:
    print("Invalid Selection")



