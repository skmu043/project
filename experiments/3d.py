import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d', adjustable='box')


x = [50,50,100]
y = [50,50,100]
z = [0,50,100]


ax.set_xlabel('X', fontsize=20, rotation=150)
ax.set_ylabel('Y')
ax.set_zlabel('Z', fontsize=30, rotation=60)


#ax.scatter(50,50,50)
#ax.scatter(51,51,51)
#ax.scatter(52,52,52)



ax.plot(x,y,z)


plt.show()
