import shelve, os
import matplotlib.pyplot as plt

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

data_dr = os.getcwd() + '/data'
data_archives = os.listdir(data_dr)

for si in data_archives:
    print(data_archives.index(si) , si)

select = int(input("Select [0-"+ str(len(data_archives) - 1)+ "]: "))

if(select <= len(os.listdir(data_dr))-1):

    s = shelve.open(data_dr + "/" + str(data_archives[select]) + "/dyke.data.db")

    try :

        time            = s['time']
        rF              = s['rF']
        rP              = s['rP']

        plot_b_p()             #plot biotic force and P

    finally:
        s.close()

else:
    print("Invalid Selection")



