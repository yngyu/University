from numpy import *
import matplotlib.pyplot as plt
import csv
from scipy.optimize import rosen, differential_evolution, minimize, basinhopping
import numpy.linalg as la

with open("satellite-temp.csv") as f:
    reader = csv.reader(f)
    l = [row for row in reader]

l = array(l).T
time = l[0]
temp = l[1]
time = [float(each_time) for each_time in time]
temp = [float(each_temp) for each_temp in temp]
time = array(time)
temp = array(temp)

diff = 1e-2

T_0 = arange(12.5,12.8,diff)
a = arange(2.8,3.0,diff)
omega = arange(1.0,1.2,diff)
theta = arange(-pi/10,pi/10,diff)

e = empty((len(T_0), len(a), len(omega), len(theta)))
"""
for i in range(len(T_0)):
    for j in range(len(a)):
        for k in range(len(omega)):
            for l in range(len(theta)):
                tempe = temp - T_0[i] - a[j]*sin(omega[k]*time + theta[l])
                tempe = tempe**2
                sume = sum(tempe)
                e[i][j][k][l] = sume

print(amin(e))
min_e = 1e3
min_i = 0
min_j = 0
min_k = 0
min_l = 0
for i in range(len(T_0)):
    for j in range(len(a)):
        for k in range(len(omega)):
            for l in range(len(theta)):
                if e[i][j][k][l] < min_e:
                    min_e = e[i][j][k][l]
                    min_i = i
                    min_j = j
                    min_k = k
                    min_l = l
print(T_0[min_i])
print(a[min_j])
print(omega[min_k])
print(theta[min_l])
"""

T_0_def = 12.67
a_def = 2.92
omega_def = 1.06
theta_def = -0.194

omega_def /= 1000
time_v = arange(0,14000,10)
sin_v = T_0_def + a_def*sin(omega_def*time_v + theta_def)

from scipy.optimize import curve_fit

def func(t,T_0,a,omega,theta):
    return T_0 + a*sin(omega*t + theta)

plt.scatter(time,temp)
plt.plot(time_v,sin_v,label = "LS")
plt.xlabel("time")
plt.ylabel("temperature")
plt.legend(loc = "lower right")
plt.savefig("460.png",bbox_inches="tight", pad_inches=0.05)
plt.show()

time /= 5000
popt, pcov = curve_fit(func,time,temp)

popt[2] /= 5000
time *= 5000
plt.scatter(time,temp)
plt.plot(time_v,sin_v,label = "LS")
plt.plot(time_v,func(time_v, *popt),label = "scipy")
print(sum((temp - func(time, *popt))**2)/len(temp))
plt.xlabel("time")
plt.ylabel("temperature")
plt.legend(loc = "lower right")
plt.savefig("461.png",bbox_inches="tight", pad_inches=0.05)
plt.show()
print(popt)


gosa_list = []
def GD(time,temp,T_0,a,omega,theta):
    step_size = 200
    n = len(time)
    step = linspace(0,step_size,step_size + 1)
    step_dif = 1e-12
    for i in range(step_size):
        dE_dT_0 = sum(2*(a*sin(omega*time + theta) + T_0 - temp))
        dE_da = sum(2*sin(omega*time + theta)*(a*sin(omega*time + theta) + T_0 - temp))
        dE_domega = sum(2*time*a*cos(omega*time + theta)*(a*sin(omega*time + theta) + T_0 - temp))
        dE_dtheta = sum(2*a*cos(omega*time + theta)*(a*sin(omega*time + theta) + T_0 - temp))
        
        T_0 -= step_dif*dE_dT_0
        a -= step_dif*dE_da
        omega -= step_dif*dE_domega
        theta -= step_dif*dE_dtheta

        gosa = sum((temp - func(time,T_0, a, omega, theta))**2)
        gosa_list.append(gosa)
    return [T_0, a, omega, theta]

init = [12.5,3,1e-3,-0.3]
popt2 = GD(time,temp,init[0],init[1],init[2],init[3])
gosa_list = array(gosa_list)
iteration_number = empty_like(gosa_list)
for i in range(len(gosa_list)):
    iteration_number[i] = i+1

"""
plt.scatter(time,temp)
#plt.plot(time_v,sin_v,label = "LS")
plt.plot(time_v,func(time_v, *popt2),label = "GD")
plt.plot(time_v,func(time_v, *popt),label = "scipy")
plt.plot(time_v,func(time_v, *init), label = "init")
plt.legend()
plt.show()
"""

gosa_list_2 = []
def Adam(time,temp,T_0,a,omega,theta):
    iteration_number = 200
    n = len(time)
    #hyper parametar
    eta = 1e-4
    rho_1 = 0.9
    rho_2 = 0.999
    epsilon = 1e-8

    m = zeros(4)
    v = zeros(4)

    w = array([T_0, a, omega,theta])

    for i in range(1,iteration_number+1):
        T_0, a, omega,theta = w

        gosa = sum((temp - func(time,T_0, a, omega, theta))**2)
        gosa_list_2.append(gosa)

        dE_dT_0 = sum(2*(a*sin(omega*time + theta) + T_0 - temp))
        dE_da = sum(2*sin(omega*time + theta)*(a*sin(omega*time + theta) + T_0 - temp))
        dE_domega = sum(2*time*a*cos(omega*time + theta)*(a*sin(omega*time + theta) + T_0 - temp))
        dE_dtheta = sum(2*a*cos(omega*time + theta)*(a*sin(omega*time + theta) + T_0 - temp))
        
        g = array([dE_dT_0, dE_da, dE_domega, dE_dtheta])
        m = rho_1*m + (1 - rho_1)*g
        v = rho_2*v + (1 - rho_2)*la.norm(g)**2
        m_hat = m/(1 - rho_1**i)
        v_hat = v/(1 - rho_2**i)
        delta_w = - eta/sqrt(v_hat + epsilon)*m_hat
        w += delta_w
    
    return w

init = [12.5,3,1e-3,-0.3] 
popt3 =  Adam(time,temp,init[0],init[1],init[2],init[3])

plt.plot(iteration_number, gosa_list, label = "GD")
plt.plot(iteration_number, gosa_list_2, label = "Adam")
plt.xlabel("iteration number")
plt.ylabel("E")
plt.legend()
plt.savefig("462.png",bbox_inches="tight", pad_inches=0.05)
plt.show()

plt.scatter(time,temp)
plt.plot(time_v,func(time_v, *popt3),label = "Adam")
plt.plot(time_v,func(time_v, *popt2),label = "GD")
plt.plot(time_v,func(time_v, *popt),label = "scipy")
plt.plot(time_v,func(time_v, *init), label = "init")
plt.xlabel("time")
plt.ylabel("temperature")
plt.legend()
plt.savefig("463.png",bbox_inches="tight", pad_inches=0.05)
plt.show()

print(sum((temp - func(time, *popt3))**2)/len(temp))