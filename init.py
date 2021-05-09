import os
structure_check = ['harness','experiments','data','analytics',]

def pre_flight_checks():
    check = True
    for s_check in structure_check:
        if(os.path.isdir(s_check)==False):
            check = False
    #print(check)

pre_flight_checks()

#print(os.listdir(os.getcwdb()))




while True:
    print (30 * "-" , "MENU" , 30 * "-")
    for s in structure_check:
        print(structure_check.index(s),s)

    print (67 * "-")
    choice = int(input("Enter your choice [0-"+ str(len(structure_check) - 1)+ "]: "))


    if(choice <= len(structure_check)-1):
        print(len(structure_check))
        #print(os.listdir(structure_check[0]))
        for s in os.listdir(structure_check[0]):
            print(os.listdir(structure_check[0]).index(s),s)
        choice = int(input("Enter your choice [0-"+ str(len(os.listdir(structure_check[0])) - 1)+ "]: "))
        if(choice <= len(os.listdir(structure_check[0]))-1):
            os.system('python simulations/dyke.py')
        else:
            print("Invalid Selection")
    else:
        print("Invalid Selection")

