import os

#"Args: K, R, P, E, start, end, step, ENumber, Niche"
#"e.g                                K=100, R=100, P=0, E=10, start=0, end=200, step=0.01, EN = 2, OE = 5"
os.system("python " + os.getcwd() + "/experiments/dyke_space.py 100 100 0 10 0 200 0.01 1 5")