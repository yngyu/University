from numpy import *
import matplotlib.pyplot as plt
import numpy.linalg as la
from scipy.optimize import rosen, differential_evolution, minimize, basinhopping

year = []
price = []

with open("mazdas_55.txt") as f:
    lines = f.readlines()
    for line in lines:
        line = line.split()
        if len(line) == 2:
            year.append(line[0])
            price.append(line[1])

year.pop(0)
price.pop(0)

year = [float(each_year) for each_year in year]
price = [float(each_price) for each_price in price]

n = len(year)
year = array(year)
price = array(price)

m = 2

def IRLS():
    X = empty((n,m+1))
    for k in range(m+1):
        X[:,k] = year**k
    y = price

    alp = la.solve(dot(X.T,X),dot(X.T,y))
    first_alp = alp

    W = zeros((n,n))

    while True:
        gosa = 1/abs(y - dot(X,alp))
        for i in range(n):
            W[i][i] = gosa[i]
        new_alp = la.solve(dot(X.T,dot(W,X)),dot(dot(X.T,W),y))
        if la.norm(new_alp - alp) < 10:
            alp = new_alp
            break
        else:
            alp = new_alp

    x_v = arange(65,95,1e-2)
    y_iv = zeros_like(x_v)
    y_ev = zeros_like(x_v)
    for k in range(len(alp)):
        y_iv += first_alp[k]*x_v**k
        y_ev += alp[k]*x_v**k

    plt.scatter(year,price)
    plt.plot(x_v,y_iv,color = "C1",label = "MSE")
    plt.plot(x_v,y_ev,color = "C2",label = "MAE")
    plt.ylim(-5000,45000)
    plt.xlabel("year")
    plt.ylabel("price")
    plt.legend(loc = "lower right")
    #plt.savefig("550.png",bbox_inches="tight", pad_inches=0.05)
    plt.savefig("552.png",bbox_inches="tight", pad_inches=0.05)
    plt.show()

    print(first_alp)
    print(alp)
    return alp

def LP():
    X = empty((n,m+1))
    for k in range(m+1):
        X[:,k] = year**k
    y = price
    first_alp = la.solve(dot(X.T,X),dot(X.T,y))

    def f(args):
        e = abs(y-dot(X,args))
        return sum(e)
    
    bounds = [(-1e9,1e9), (-1e9,1e9), (-1e9,1e9)]

    ret = differential_evolution(f, bounds)
    alp = ret["x"]

    x_v = arange(65,95,1e-2)
    y_iv = zeros_like(x_v)
    y_ev = zeros_like(x_v)
    for k in range(len(alp)):
        y_iv += first_alp[k]*x_v**k
        y_ev += alp[k]*x_v**k

    plt.scatter(year,price)
    plt.plot(x_v,y_iv,color = "C1")
    plt.plot(x_v,y_ev,color = "C2")
    plt.ylim(-5000,45000)
    plt.show()

    print(first_alp)
    print(alp)

    return alp

X = empty((n,m+1))
for k in range(m+1):
    X[:,k] = year**k
y = price
alp = la.solve(dot(X.T,X),dot(X.T,y))

IRLS_alp = IRLS()
LP_alp = LP()

x_v = arange(65,95,1e-2)
y_iv = zeros_like(x_v)
y_IRLS = zeros_like(x_v)
y_LP = zeros_like(x_v)
for k in range(len(alp)):
    y_iv += alp[k]*x_v**k
    y_IRLS += IRLS_alp[k]*x_v**k
    y_LP += LP_alp[k]*x_v**k

plt.scatter(year,price)
plt.plot(x_v,y_iv,color = "C1",label = "MSE")
plt.plot(x_v,y_IRLS,color = "C2",label = "IRLS")
plt.plot(x_v,y_LP,color = "C3",label = "differential evolution")
plt.xlabel("year")
plt.ylabel("price")
plt.legend()
#plt.savefig("551.png",bbox_inches="tight", pad_inches=0.05)
plt.savefig("553.png",bbox_inches="tight", pad_inches=0.05)
plt.ylim(-5000,45000)
plt.show()