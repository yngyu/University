def f(a,b,c):
    return b
def g(a,b,c):
    return c
def h(a,b,c):
    return (-1)*a*c/2

sup = 1
inf = 0

while True:
    x = [0]
    y = [0]
    z = [(sup+inf)/2]

    delta_t = 0.0001
    cal_number = 100000

    for i in range(1,cal_number):
        k0 = delta_t*f(x[-1],y[-1],z[-1])
        l0 = delta_t*g(x[-1],y[-1],z[-1])
        m0 = delta_t*h(x[-1],y[-1],z[-1])
        k1 = delta_t*f(x[-1]+k0/2,y[-1]+l0/2,z[-1]+m0/2)
        l1 = delta_t*g(x[-1]+k0/2,y[-1]+l0/2,z[-1]+m0/2)
        m1 = delta_t*h(x[-1]+k0/2,y[-1]+l0/2,z[-1]+m0/2)
        k2 = delta_t*f(x[-1]+k1/2,y[-1]+l1/2,z[-1]+m1/2)
        l2 = delta_t*g(x[-1]+k1/2,y[-1]+l1/2,z[-1]+m1/2)
        m2 = delta_t*h(x[-1]+k1/2,y[-1]+l1/2,z[-1]+m1/2)
        k3 = delta_t*f(x[-1]+k2,y[-1]+l2,z[-1]+m2)
        l3 = delta_t*g(x[-1]+k2,y[-1]+l2,z[-1]+m2)
        m3 = delta_t*h(x[-1]+k2,y[-1]+l2,z[-1]+m2)
        x.append(x[-1]+(k0+2*k1+2*k2+k3)/6)
        y.append(y[-1]+(l0+2*l1+2*l2+l3)/6)
        z.append(z[-1]+(m0+2*m1+2*m2+m3)/6)
    
    if 0.999 <= y[cal_number-1] and y[cal_number-1] <= 1.001:
        break
    elif y[cal_number-1] > 1:
        sup = (sup+inf)/2
    elif y[cal_number-1] < 1:
        inf = (sup+inf)/2

print("f''(0)")
print((sup+inf)/2)
print("f'(10)")
print(y[cal_number-1])

import numpy as np
import matplotlib.pyplot as plt
import math

t = []
for i in range(0,cal_number):
    t.append(i*delta_t)
p = np.array(t)
q = np.array(y)
plt.plot(q,p)
plt.xlabel("f'(η)")
plt.ylabel("η")
plt.show()