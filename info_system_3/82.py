import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

I_x = 1.9
I_y = 1.6
I_z = 2.0
I = np.array([I_x, I_y, I_z])
diff = 1e-2
sigma_w = 1e-2
sigma_v = 4e-2

def Euler(u,noise_flag):
    x,y,z = u
    dxdt = (I_y - I_z)/I_x*y*z
    dydt = (I_z - I_x)/I_y*x*z
    dzdt = (I_x - I_y)/I_z*x*y
    if noise_flag:
        w = np.random.normal(loc = 0.0,scale = sigma_w, size = 3)
    else:
        w = np.zeros(3)
    return u + diff*(np.array([dxdt, dydt, dzdt]) + w/I)

u_init = np.array([0.1, 17*2*np.pi/60 + 0.1, 0])
u = [u_init]

t = np.arange(0,40,diff)

for i in range(len(t)-1):
    u.append(Euler(u[-1],True))

u = np.array(u)

v = np.random.normal(loc = 0.0, scale = sigma_v, size = (u.shape[0], u.shape[1]))

y = u + v

def Kalman_Filter(y_oberve):
    y_kalman = [y_oberve[0]]
    m = y_oberve[0]
    temp = np.full(3,sigma_w**2)
    V = np.diag(temp) #予測値の分散行列、テキトーに置いた
    temp = np.array([sigma_w**2/I_x**2, sigma_w**2/I_y**2, sigma_w**2/I_z**2]) #xにかかるノイズの分散行列
    Q = np.diag(temp)

    Ones = np.ones((3,3))
    Identity = np.identity(3)
    C = Identity #観測の行列
    temp = np.full(3,sigma_v**2)
    R = np.diag(temp) #観測ノイズの分散行列
    
    for i in range(1,u.shape[0]):
        #予測
        A = np.zeros((3,3))
        x,y,z = y_oberve[i]
        A[0][1] = (I_y - I_z)/I_x*z
        A[0][2] = (I_y - I_z)/I_x*y
        A[1][0] = (I_z - I_x)/I_y*z
        A[1][2] = (I_z - I_x)/I_y*x
        A[2][0] = (I_x - I_y)/I_z*y
        A[2][1] = (I_x - I_y)/I_z*x
        A *= diff
        A += Ones

        m = Euler(m,False)
        V = np.dot(A,np.dot(V,A.T)) + Q

        #修正
        P = la.inv(np.dot(C,np.dot(V, C.T) + R))
        K = np.dot(V,np.dot(C.T,P)) #カルマンゲイン
        m += np.dot(K,y_oberve[i] - np.dot(C,m))
        y_kalman.append(m)
        V = np.dot((Identity - C),V)
    
    y_kalman = np.array(y_kalman)
    return y_kalman

y_kalman = Kalman_Filter(y)

for i in range(u.shape[1]):
    plt.plot(t,y[:,i], label = "observe")
    plt.plot(t,y_kalman[:,i], label = "Kalman_Filter")
    plt.plot(t,u[:,i], label = "True")
    plt.xlabel("t[s]")
    if i == 0:
        plt.ylabel(r"$\omega_x$")
    elif i == 1:
        plt.ylabel(r"$\omega_y$")
    else:
        plt.ylabel(r"$\omega_z$")
    plt.legend(loc = "lower right")
    if i == 0:
        plt.savefig("82x.png",bbox_inches="tight", pad_inches=0.05)
    elif i == 1:
        plt.savefig("82y.png",bbox_inches="tight", pad_inches=0.05)
    else:
        plt.savefig("82z.png",bbox_inches="tight", pad_inches=0.05)
    plt.show()