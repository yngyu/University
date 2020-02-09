from scipy.stats import poisson
import matplotlib.pyplot as plt
from numpy import *

def problem1():
    mu = 30/1000

    number = arange(0,11,1)

    poisson_result = poisson.pmf(number,mu)

    plt.plot(number,poisson_result)
    plt.show()

    poisson_sum = empty_like(poisson_result)
    poisson_sum[0] = poisson_result[0]
    for i in range(1,len(poisson_sum)):
        poisson_sum[i] = poisson_sum[i-1] + poisson_result[i]

    plt.plot(number,poisson_sum)
    plt.show()

    print(poisson_result[0])

def l_func(mu):
    p1 = mu*exp(-mu)
    p2 = mu**2*exp(-mu)/2
    p3 = mu**3*exp(-mu)/6
    return pow(p1,5)*pow(p2,4)*p3

def ll_func(mu):
    p1 = mu*exp(-mu)
    p2 = mu**2*exp(-mu)/2
    p3 = mu**3*exp(-mu)/6
    p1 = log(p1)
    p2 = log(p2)
    p3 = log(p3)
    return 5*p1 + 4*p2 + p3

def problem3():
    diff = 1e-4
    mu = arange(0,3,diff)
    plt.plot(mu,l_func(mu))
    plt.scatter(1.6,l_func(1.6))
    plt.show()
    plt.plot(mu,ll_func(mu))
    plt.scatter(1.6,ll_func(1.6))
    plt.show()

problem3()