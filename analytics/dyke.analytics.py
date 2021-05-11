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

        #plot_alphas()          #plot abundance of species over temperature
        #plot_w()               #plot affects values for each species
        #plot_u()               #plot ideal growing temperature for each species
        #plot_aot()             #plot abundance of each species over time
        #plot_aot_scaled()      #plot abundance of each species over time scaled by R
        plot_aot_inc_dec()     #plot species that increase temperature and decrease temperature
        #plot_b_p()             #plot biotic force and P
        #plot_e()               #plot temperature value over time
        #plot_efp()             #plot temperature, biotic force and P over time
    finally:
        s.close()

else:
    print("Invalid Selection")



