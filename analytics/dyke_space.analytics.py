import shelve, os
import matplotlib.pyplot as plt
import sys
print(sys.version)

#plot affects values for each species
def plot_w():

    plt.figure(figsize=(20,5))
    plt.ylim(-1, 1)
    plt.title('Affects values for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Affects', fontsize=18)

    for list in w:
        plt.plot(list, 'k.', label='w')
    plt.legend(loc=5, prop={'size': 30})
    plt.show()

#plot ideal growing temperature for each species
def plot_u():
    plt.figure(figsize=(20,5))
    plt.ylim(0, R)
    plt.title('Ideal Growing Temperature for each species', fontsize=20)
    plt.xlabel('Species', fontsize=18)
    plt.ylabel('Temperature', fontsize=18)
    for list in u:
        plt.plot(list, 'k.', label='u')
    plt.legend(loc=5, prop={'size': 30})
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
    for list in rE:
        plt.plot(time,list, 'r.', label='E')
    plt.show()

#plot temperature value over time
def plot_e():
    plt.figure(figsize=(20,10))
    plt.title('Environment Variable E', fontsize=20)
    plt.xlabel('Time', fontsize=20)
    plt.ylabel('Temperature', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    for list in rE:
        plt.plot(time,list, 'r-', label='E')
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
    #for list in rF:
    #    plt.plot(time,list, 'g.', label = 'biotic force')
    for list in rE:
        plt.plot(time,list, 'r.',label = 'temperature')
    plt.legend(loc='lower right', prop={'size': 30})
    plt.axhline(y=R)

    xtrajectory=[]
    ytrajectory=[]
    for idx in range(len(rE[0])):

        xtrajectory.append(rE[0][idx])
        ytrajectory.append(rE[1][idx])

        if (idx > 0) and (idx % (end/step) == 0):
            #print(idx)
            #print(xtrajectory)
            #print(ytrajectory)
            plt.plot(xtrajectory, '-',ytrajectory ,'.', label='traj')
            del xtrajectory[:]
            del ytrajectory[:]

    plt.show()

def plot_stable():
    #plt.figure(figsize=(30,20))
    plt.title('Stable Points', fontsize=40)
    #plt.xlabel('E Values', fontsize=40)
    #plt.ylabel('E Values', fontsize=40)
    #plt.xticks(fontsize=20)
    #plt.yticks(fontsize=20)
    #plt.ylim(-50, R+20)
    #plt.xlim(0, end)

    #print(rE)

    xtrajectory=[]
    ytrajectory=[]
    ztime = []

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(0,2000)
    ax.set_ylim3d(-50,100)
    ax.set_zlim3d(0,100)

    index_z = 0
    for idx in range(len(rE[0])):

        xtrajectory.append(rE[0][idx])
        ytrajectory.append(rE[1][idx])
        ztime.append(index_z)
        index_z += 1

        if (idx > 0) and (idx % (end/step) == 0):
            #print(idx)
            #print(xtrajectory)
            # print(ytrajectory)
            #print(ztime)
            #plt.plot(xtrajectory,ytrajectory, '-', label='traj')
            if(rE[0][idx]<=100 and rE[1][idx]>=0):
                print("Stable Point : Start Temp: ", xtrajectory[0], ytrajectory[0])
            ax.plot(ztime, xtrajectory, ytrajectory)
            xtrajectory.clear()
            ytrajectory.clear()
            ztime.clear()
            index_z = 0

    plt.legend(loc='lower right', prop={'size': 30})
    plt.show()


data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):
    #print(os.listdir(data_dr))
    #print(data_dr + "/" + str(os.listdir(data_dr)[select]) + "/dyke.data.db")
    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke_space.data")
    #print(data_dr + "/" + str(data_archives[select]) + "/dyke_space.data.db")
    #print("/Users/sumeet.kumar/IdeaProjects/project/data/1624411616.2321012.dyke_space/dyke_space.data")
    #s = shelve.open("/Users/sumeet.kumar/IdeaProjects/project/data/1624411616.2321012.dyke_space/dyke_space.data")
    print(s)
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
        rAx             = s['rAx']

        K               = s['K']
        R               = s['R']
        E               = s['E']
        start           = s['start']
        end             = s['end']
        step            = s['step']
        N               = s['N']
        OEn             = s['OEn']



        print("0 plot affects values for each species")
        print("1 plot ideal growing temperature for each species")
        print("2 plot abundance of each species over time")
        print("3 plot abundance of each species over time scaled by R")
        print("4 plot temperature value over time")
        print("5 plot temperature, biotic force")
        print("7 Stable Points")
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
            plot_e()               #plot temperature value over time
        elif select == 5:
            plot_efp()             #plot temperature, biotic force over time
        elif select == 7:
            plot_stable()         # [E1 values over time >>>>>]
                                  # [E2 values over time >>>>>]
                                  #.
                                  #.
                                  # [En values over time >>>>>]
        else:
            print("Invalid Selection")
    finally:
        s.close()

else:
    print("Invalid Selection")
