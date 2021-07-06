# Python 3.9.1
import math, random, sys, os, shelve, time
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d', adjustable='box')


x = [50,50,100]
y = [50,50,100]
z = [0,50,100]

ax.set_xlabel('X', fontsize=20, rotation=150)
ax.set_ylabel('Y')
ax.set_zlabel('Z', fontsize=30, rotation=60)
ax.plot(x,y,z)
plt.show()


#ax.scatter(50,50,50)
#ax.scatter(51,51,51)
#ax.scatter(52,52,52)

E       = []
alpha   = []
start   = 0
end     = 100
step    = 1

check = ("","")
checkone = []

for xtime in np.arange (start, end, step):
    E.append(xtime)
    alpha.append(round((math.e) ** ((-1) * (((abs((E[-1])-50)) ** 2) / (2*(5**2)))), 5))

    check = (E[-1],alpha[-1])
    checkone.append(check)
#print(E)
#print(alpha)

print(checkone)

plt.figure()
plt.ylim(-1, 1)
plt.title('Abundance Check', fontsize=20)
plt.xlabel('Temp', fontsize=18)
plt.ylabel('Alpha', fontsize=18)
plt.plot(E, alpha)
plt.legend(loc=5, prop={'size': 30})
plt.show()



