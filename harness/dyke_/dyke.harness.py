import os, sys

#"Args: K, R, P, E, start, end, step"
#"e.g K=100, R=100, P=0, E=10, start=0, end=200, step=0.01"
print("Python version: ", sys.version, "Version info: ", sys.version_info)


os.system("python3.9 " + os.getcwd() + "/experiments/dyke.py 100 100 0 10 0 200 0.01")

