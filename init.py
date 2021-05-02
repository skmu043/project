import os

structure_check = {'analytics','harness','results_data','simulations'}

def pre_flight_checks():
    check = True
    for s_check in structure_check:
        if(os.path.isdir(s_check)==False):
            check = False
    print(check)

pre_flight_checks()



print(os.listdir(os.getcwdb()))
