from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt


def fun(x):
    return [x[0]  + 0.5 * (x[0] - x[1])**3 - 1.0,0.5 * (x[1] - x[0])**3 + x[1]]

def jac(x):
    return np.array([[1 + 1.5 * (x[0] - x[1])**2,-1.5 * (x[0] - x[1])**2],[-1.5 * (x[1] - x[0])**2,1 + 1.5 * (x[1] - x[0])**2]])

sol = optimize.root(fun, [0, 0], jac=jac, method='hybr')

print(sol.x)


def f1(x):
    #return((x**3) + (2*(x**2)) - (2*x) - 5)
    return(x**2 -1000)

x = []
y = []

X1 = -100
Y1 = 100

zeroz = []

for xi in np.arange(X1, Y1, 0.001):
    x.append(xi)
    y.append(f1(xi))
    #print(y[-1])
    if(y[-1] == 0):
        zeroz.append(x[-1])

sol = optimize.root(f1, [-100, 100], jac=False, method='hybr')

print(sol.x)

plt.figure(figsize=(20,10))
plt.title('xy', fontsize=40)
plt.xlabel('x', fontsize=40)
plt.ylabel('y', fontsize=40)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.ylim(X1, Y1)
#plt.xlim(0, end)
for x3 in sol.x:
    print(x3)
    plt.axvline(x=x3)

plt.axvline(x=0)
plt.axhline(y=0)

plt.plot(x,y, 'r--',label = 'roots')
plt.show()


