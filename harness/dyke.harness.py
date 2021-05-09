import os

print(os.getcwd())

os.system("python " + os.getcwd() + "/experiments/dyke.py 100 100 0 10 0 200 0.01")