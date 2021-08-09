import numpy as np
import matplotlib.pyplot as plt

# Creating data set
w = 3
Y, X = np.mgrid[-w:w:100j, -w:w:100j]
U = -1 - X**2 + Y
V = 1 + X - Y**2
speed = np.sqrt(U**2 + V**2)

# Creating plot
fig = plt.figure(figsize = (12, 7))
#plt.streamplot(X, Y, U, V, density = 1)

# show plot
#plt.show()


# Creating dataset
x = np.arange(0, 2)
y = np.arange(0, 2)

# Creating grids
X, Y = np.meshgrid(x, y)

# x-component to the right
u = np.ones((2, 2))

# y-component zero
v = np.zeros((2, 2))

fig = plt.figure(figsize = (12, 7))

print(x)
print(y)
print(X)
print(Y)
print(u)
print(v)

# Plotting stream plot
plt.streamplot(X, Y, u, v, density = 0.5)

# show plot
plt.show()