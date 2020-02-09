from numpy import *
import matplotlib.pyplot as plt
import numpy.linalg as la
from scipy.optimize import rosen, differential_evolution, minimize, basinhopping
from scipy.stats import norm
import csv

with open("ballthrow.csv") as f:
    reader = csv.reader(f)
    l = []
    for row in reader:
        l.append(float(row[0]))

l = array(l)
l_min = min(l)
l_max = max(l)

n = len(l)
data = l
def dif_evo():
    def func(args):
        alpha_M, mu_M, sigma_M, mu_F, sigma_F = args
        man = alpha_M/sqrt(2*pi*sigma_M)*exp(-(l - mu_M)**2/(2*sigma_M))
        woman = (1-alpha_M)/sqrt(2*pi*sigma_F)*exp(-(l - mu_F)**2/(2*sigma_F))
        temp = sum(log(man + woman))
        return -temp

    bounds = [(0.05,0.95), (25,50), (1,1000),(10,30), (1,1000)]

    ret = differential_evolution(func, bounds)

    print(ret)

    param = ret["x"]
    alpha_M, mu_M, sigma_M, mu_F, sigma_F = param
    diff = 1e-3
    x_v = arange(0,60,diff)
    a = 5
    man_v = a*n*alpha_M/sqrt(2*pi*sigma_M)*exp(-(x_v - mu_M)**2/(2*sigma_M))
    woman_v = a*n*(1-alpha_M)/sqrt(2*pi*sigma_F)*exp(-(x_v - mu_F)**2/(2*sigma_F))

    plt.plot(x_v,man_v,label = "male")
    plt.plot(x_v,woman_v,label = "female")
    plt.plot(x_v,man_v + woman_v,label = "all")
    plt.hist(l, bins = 16,color = "gray")
    plt.savefig("470.png",bbox_inches="tight", pad_inches=0.05)
    plt.legend()
    plt.show()

def EM_algorithm():
    #alpha_M, mu_M, mu_F, sigma_M, sigma_F
    init_list = [0.65,35,10,100,10]
    args = init_list

    iteration_number = 10000

    def calculate_w(j,args): #iはどのdataか,jはmaleかfemaleか
        alpha_M, mu_M, mu_F, sigma_M, sigma_F = args
        temp_male = alpha_M*norm.pdf(x = data, loc = mu_M, scale = sqrt(sigma_M))
        temp_female = (1-alpha_M)*norm.pdf(x = data, loc = mu_F, scale = sqrt(sigma_F))

        if j == 0:
            return temp_male/(temp_male + temp_female)
        else:
            return temp_female/(temp_male + temp_female)
    
    def calculate_large_N_K(j,args):
        return sum(calculate_w(j,args))
    
    def calculate_new_mu(j,args):
        temp = sum(calculate_w(j,args)*data)
        temp /= calculate_large_N_K(j,args)
        return temp
    
    def calculate_new_sigma2(j,args):
        if j == 0:
            mu = args[1]
        else:
            mu = args[2]
        temp = sum(calculate_w(j,args)*(data - mu)**2)
        temp /= calculate_large_N_K(j,args)
        return temp

    def calculate_new_alpha_M(args):
        return  calculate_large_N_K(0,args)/n

    def calculate_LL(args):
        alpha_M, mu_M, mu_F, sigma_M, sigma_F = args
        temp_male = alpha_M*norm.pdf(x = data, loc = mu_M, scale = sqrt(sigma_M))
        temp_female = (1-alpha_M)*norm.pdf(x = data, loc = mu_F, scale = sqrt(sigma_F))

        return sum(log(temp_male + temp_female))

    def iteration(old_args):
        LL_list = []
        args = old_args
        LL_list.append(calculate_LL(args))
        for i in range(iteration_number):
            new_alpha_M = calculate_new_alpha_M(args)
            new_mu_M = calculate_new_mu(0,args) #0がmale
            new_mu_F = calculate_new_mu(1,args) #1がfemale
            new_sigma_M = calculate_new_sigma2(0,args)
            new_sigma_F = calculate_new_sigma2(1,args)

            new_args = [new_alpha_M, new_mu_M, new_mu_F, new_sigma_M, new_sigma_F]

            LL_list.append(calculate_LL(new_args))

            if la.norm(array(args) - array(new_args)) < 1e-6:
                args = new_args
                print(i)
                break
            else:
                args = new_args
                
        return args,LL_list
    
    alp,LL = iteration(args) 
    print(alp)
    LL_number = empty(len(LL))
    for i in range(len(LL)):
        LL_number[i] = i+1
    plt.plot(LL_number,LL)
    plt.xlabel("iteration")
    plt.ylabel("LL")
    plt.savefig("471.png",bbox_inches="tight", pad_inches=0.05)
    plt.show()

    alpha_M, mu_M, mu_F, sigma_M, sigma_F = alp
    diff = 1e-3
    x_v = arange(0,60,diff)
    a = 5

    man_v = a*n*alpha_M*norm.pdf(x = x_v, loc = mu_M, scale = sqrt(sigma_M))
    woman_v = a*n*(1-alpha_M)*norm.pdf(x = x_v, loc = mu_F, scale = sqrt(sigma_F))

    plt.plot(x_v,man_v,label = "male")
    plt.plot(x_v,woman_v,label = "female")
    plt.plot(x_v,man_v + woman_v,label = "all")
    plt.hist(l, bins = 16,color = "gray")
    plt.savefig("472.png",bbox_inches="tight", pad_inches=0.05)
    plt.legend()
    plt.show()

dif_evo()
EM_algorithm()