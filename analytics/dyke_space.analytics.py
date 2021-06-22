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

#plot ideal growing temperature for each species
def plot_u():
    plt.figure(figsize=(20,5))
    plt.ylim(0, R)
    plt.title('Ideal Growing Temperature for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    plt.plot(u, 'k.', label='u')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

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

#plot abundance of each species over time
def plot_aot():
    plt.figure(figsize=(20,10))
    plt.title('Abundance over Time', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Abundance', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for x in range(K):
        plt.plot(time,rAx[x],label = 'id %s'%x)
    plt.show()

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

#plot biotic force and P
def plot_b_p():
    plt.figure(figsize=(20,10))
    plt.title('Biotic Force and P', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Value for Biotic Force and P', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot(time,rF, 'g-', label='F')
    plt.plot(time,rP, 'b--', label='P')
    plt.show()

#plot temperature value over time
def plot_e():
    plt.figure(figsize=(20,10))
    plt.title('Environment Variable E', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.plot(time,rE, 'r-', label='E')
    plt.axhline(y=R)
    plt.show()

#plot temperature, biotic force and P over time
def plot_efp():
    plt.figure(figsize=(30,20))
    plt.title('Simulation Values over Time', fontsize=40)
    plt.xlabel('Time', fontsize=40)
    plt.ylabel('Values', fontsize=40)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylim(-50, R+20)
    plt.xlim(0, end)
    plt.plot(time,rF, 'g-', label = 'biotic force')
    plt.plot(time,rP, 'b--', label = 'perturbing force(rate)')
    plt.plot(time,rE, 'r-',label = 'temperature')
    plt.plot(time,rEt, 'k.',label = 'temperature (without biota)')
    plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)
    plt.show()


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):
    #print(os.listdir(data_dr))
    #print(data_dr + "/" + str(os.listdir(data_dr)[select]) + "/dyke.data.db")
    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke_space.data.db")

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



        print("0 plot affects values for each species")
        print("1 plot ideal growing temperature for each species")
        print("2 plot abundance of each species over time")
        print("3 plot abundance of each species over time scaled by R")
        print("4 plot species that increase temperature and decrease temperature")
        print("5 plot biotic force and P")
        print("6 plot temperature value over time")
        print("7 plot temperature, biotic force and P over time")

        select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

        if select == 0:
            plot_w()               #plot affects values for each species
        elif select == 1:
            plot_u()               #plot ideal growing temperature for each species
        elif select == 2:
            plot_aot()             #plot abundance of each species over time
        elif select == 3:
            plot_aot_scaled()      #plot abundance of each species over time scaled by R
        elif select == 4:
            plot_aot_inc_dec()      #plot species that increase temperature and decrease temperature
        elif select == 5:
            plot_b_p()             #plot biotic force and P
        elif select == 6:
            plot_e()               #plot temperature value over time
        elif select == 7:
            plot_efp()             #plot temperature, biotic force and P over time
        else:
            print("Invalid Selection")
    finally:
        s.close()

else:
    print("Invalid Selection")
