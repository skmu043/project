import shelve, os
import matplotlib.pyplot as plt

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

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):

    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke.data.db")

    try :

        time            = s['time']
        rE              = s['rE']
        R               = s['R']

        plot_e()               #plot temperature value over time

    finally:
        s.close()

else:
    print("Invalid Selection")



