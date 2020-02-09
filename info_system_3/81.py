from numpy import *
import numpy.random as ra
import matplotlib.pyplot as plt
import seaborn as sns

mu = [20,10]
mu = array(mu)
sigma = [[30,20],[20,50]]

n = 10000

values = ra.multivariate_normal(mu,sigma,n)

A = [[sqrt(sigma[0][0]),0],[sigma[0][1]/sqrt(sigma[0][0]),sqrt(sigma[1][1] - sigma[0][1]**2/sigma[0][0])]]
A = array(A)

values_2 = [ra.randn(n),ra.randn(n)]
values_2 = array(values_2)

temp = dot(A,values_2)

values_2[0] = temp[0] + mu[0]
values_2[1] = temp[1] + mu[1]

sns.jointplot(values[:,0],values[:,1])
plt.savefig("810.png",bbox_inches="tight", pad_inches=0.05)
sns.jointplot(values_2[0],values_2[1])
plt.savefig("811.png",bbox_inches="tight", pad_inches=0.05)
plt.show()