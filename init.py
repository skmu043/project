import os
structure_check = ['analytics','harness','results_data','simulations']

def pre_flight_checks():
    check = True
    for s_check in structure_check:
        if(os.path.isdir(s_check)==False):
            check = False
    print(check)

pre_flight_checks()


print(os.listdir(os.getcwdb()))

def print_menu(): 
    print (30 * "-" , "MENU" , 30 * "-")
   
    for s in structure_check:
        print(structure_check.index(s),s)

    print (67 * "-")
  

while True:          
    print_menu()   
    choice = int(input("Enter your choice [1-4]: "))
    print(choice)

    if(choice <= len(structure_check)-1):
        print(len(structure_check))
        print(os.listdir(structure_check[0]))
    else:
        print("Invalid Selection")
        os.system('python simulations/dyke.py')
