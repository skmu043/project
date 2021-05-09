import os
structure_ = ['harness','experiments','data','analytics',]

def pre_flight_checks():
    check = True
    for si in structure_:
        if(os.path.isdir(si)==False):
            check = False

pre_flight_checks()
#print(os.listdir(os.getcwdb()))

while True:
    for si in structure_:
        print(structure_.index(si),si)

    select = int(input("Select [0-"+ str(len(structure_) - 1)+ "]: "))

    if(select <= len(structure_)-1):
        #print(len(structure_))
        #print(os.listdir(structure_check[0]))
        for si in os.listdir(structure_[select]):
            print(os.listdir(structure_[select]).index(si),si)
        run_select = int(input("Select [0-"+ str(len(os.listdir(structure_[select])) - 1)+ "]: "))
        if(run_select <= len(os.listdir(structure_[select]))-1):
            run_= "python " + structure_[select] + "/" +(os.listdir(structure_[select]))[run_select]
            print(run_)
            os.system(run_)
        else:
            print("Invalid Selection")
    else:
        print("Invalid Selection")

