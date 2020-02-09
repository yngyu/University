from numpy import *
import matplotlib.pyplot as plt
import numpy.linalg as la
import csv

distance = []
velocity = []

with open("Hubbles_constant.csv") as f:
    reader = csv.reader(f)
    for row in reader:
        distance.append(row[0])
        velocity.append(row[1])

distance.pop(0)
velocity.pop(0)

distance = [float(each) for each in distance]
velocity = [float(each) for each in velocity]

distance = array(distance)
velocity = array(velocity)

a_1 = sum(distance*velocity)/sum(distance**2)

x_v = arange(-0.5,2.5,1e-4)
y_v = a_1*x_v

X_m = empty((len(distance),2))
X_m[:,0] = 1
X_m[:,1] = distance

Y = velocity

alp = la.solve(dot(X_m.T,X_m), dot(X_m.T, Y))

y_a = alp[0] + alp[1]*x_v
""""
plt.scatter(distance, velocity)
plt.plot(x_v,y_v, label = r"$y = ax$")
plt.plot(x_v,y_a, label = r"$y = ax + b$")
plt.xlabel("distance")
plt.ylabel("velocity")
plt.legend(loc = "lower right")
plt.savefig("53.png",bbox_inches="tight", pad_inches=0.05)
plt.show()
"""
AIC = empty(2)

n = len(distance)
print(n)

AIC[0] = sum((velocity - a_1*distance)**2)
AIC[1] = sum((velocity - (alp[0] + alp[1]*distance))**2)

AIC[0] = n*log(AIC[0]/n) + 2*n/(n-2-1)
AIC[1] = n*log(AIC[1]/n) + 2*2*n/(n-2*2-1)
print(a_1)
print(alp)
print(AIC)