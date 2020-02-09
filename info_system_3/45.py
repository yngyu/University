from numpy import *
import matplotlib.pyplot as plt
import numpy.linalg as la
from scipy.stats import norm
from mpl_toolkits.mplot3d import Axes3D

x = [18.2,25.82,18.24,28.6,31.1,33.6,40.46,28.27,20.1,27.91,26.18,22.12]
y = [17.05,19.8,15.98,22.07,22.83,24.55,27.27,23.57,13.58,22.8,20.3,16.59]

X_m = empty((len(x),2))
for i in range(len(x)):
    X_m[i][0] = 1
    X_m[i][1] = x[i]

y_v = array(y).T

ans = la.solve(dot(X_m.T,X_m),dot(X_m.T,y_v))

diff = 1e-2

x_a = arange(0,50,diff)
y_a = ans[0] + ans[1]*x_a

print(ans)

plt.scatter(x,y)
plt.plot(x_a,y_a, label = "LS")
plt.xlim(0,50)
plt.ylim(0,30)
plt.xlabel("Smoking rate")
plt.ylabel("Death rate")
plt.legend()
plt.show()

a = arange(0,1,diff)
b = arange(4,6,diff)
aa,bb = meshgrid(a,b)

arg_max_aa = 0
arg_max_bb = 0

def plot_l_func():
    def l_func(a,b):
        x_lv = array(x)
        y_lv = array(y)
        e = y_lv - (a*x_lv+b)
        e = e**2
        p = 1/sqrt(2*pi)*exp(-e/2)
        cal_p = prod(p)
        return cal_p

    l = empty_like(aa)
    for i in range(aa.shape[0]):
        for j in range(aa.shape[1]):
            l[i][j] = l_func(aa[i][j],bb[i][j])
    print(l.max())
    i = argmax(l)//l.shape[1]
    j = argmax(l)%l.shape[1]
    print(i,j)
    print(aa[i][j],bb[i][j])
    global arg_max_aa
    arg_max_aa = aa[i][j]
    global arg_max_bb
    arg_max_bb = bb[i][j]
    
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(aa,bb,l)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("L")
    plt.show()

def plot_ll_func():
    def ll_func(a,b):
        x_lv = array(x)
        y_lv = array(y)
        e = y_lv - (a*x_lv+b)
        e = e**2
        sum_e = sum(e)
        p = -len(e)/2*(log(2*pi)) - 1/2*sum_e
        return p

    ll = empty_like(aa)
    for i in range(aa.shape[0]):
        for j in range(aa.shape[1]):
            ll[i][j] = ll_func(aa[i][j],bb[i][j])
    
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel("a")
    ax.set_ylabel("b")
    ax.set_zlabel("L")
    ax.plot_surface(bb,aa,ll)
    plt.show()
    print(ll.max())
    i = argmax(ll)//ll.shape[1]
    j = argmax(ll)%ll.shape[1]
    print(i,j)
    print(aa[i][j],bb[i][j])
    global arg_max_aa
    arg_max_aa = aa[i][j]
    global arg_max_bb
    arg_max_bb = bb[i][j]

plot_l_func()
plot_ll_func()
y_a2 = x_a*arg_max_aa + arg_max_bb

plt.scatter(x,y,color = "C0")
plt.plot(x_a,y_a, label = "LS", color = "C1")
plt.plot(x_a,y_a2, label = "MLE", color = "C2")
plt.xlabel("Smoking rate")
plt.ylabel("Death rate")
plt.xlim(0,50)
plt.ylim(0,30)
plt.legend()
plt.show()

y_g = sum((array(y) - (ans[0] + ans[1]*array(x)))**2)
y_g2 = sum((array(y) - (arg_max_bb + arg_max_aa*array(x)))**2)
print("e_ls", y_g)
print("e_mle", y_g2)

sigma2 = arange(diff,5,1e-6)

def plot_sigma_l_func(a,b,sigma):
    x_lv = array(x)
    y_lv = array(y)
    e = y_lv - (a*x_lv+b)
    e = e**2
    p = empty_like(e)
    cal_p = ones_like(sigma)
    for j  in range(len(sigma)):
        for i in range(len(p)):
            p[i] = 1/sqrt(2*pi*sigma[j])*exp(-e[i]/2/sigma[j])
        for i in range(len(p)):
            cal_p[j] *= p[i]
    plt.plot(sigma,cal_p)
    plt.xlabel(r"$\sigma^2$")
    plt.ylabel("L")
    plt.show()
    sigma2 = sigma[argmax(cal_p)]
    print(sigma2)
    print(sqrt(sigma2))

def plot_sigma_ll_func(a,b,sigma):
    x_lv = array(x)
    y_lv = array(y)
    e = y_lv - (a*x_lv+b)
    e = e**2
    sum_e = sum(e)
    p = -len(e)/2*(log(2*pi) + log(sigma)) - 1/2/sigma*sum_e
    plt.plot(sigma,p)
    plt.xlabel(r"$\sigma^2$")
    plt.ylabel("LL")
    plt.show()
    sigma2 = sigma[argmax(p)]
    print(sigma2)
    print(sqrt(sigma2))

plot_sigma_ll_func(ans[1],ans[0],sigma2)