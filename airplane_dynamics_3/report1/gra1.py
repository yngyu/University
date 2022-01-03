import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m = 7500
S = 27.9
g = 9.8
rho = 0.736423
U = 180

C_lalpha = 4.30
C_d0 = 0.0548
K = 3.02
alpha = 0.050617000000425585
L = rho*U**2*S/2*C_lalpha*alpha
D = rho*U**2*S/2*(C_d0+K*alpha**2)
T = (20841.80471039165+20842.199605483533)/2

Phi_0 = 30/180*np.pi
omega = 2*np.pi/40

dif = 1e-3
t = np.arange(0,80+dif,dif)
Phi = np.empty_like(t)

for i in range(0,t.size):
    if t[i] <= 20:
        Phi[i] = 0
    elif t[i] > 20 and t[i] <= 60:
        Phi[i] = Phi_0*np.sin(omega*(t[i]-20))
    else:
        Phi[i] = 0

x_E = np.empty((6,t.size))
U_mat = np.empty((2,t.size))
angle = np.empty((4,t.size))
x_E.T[0] = 0
x_E[5][0] = 5000
U_mat[0][0] = 0
U_mat[1][0] = U

angle.T[0] = 0

for i in range(1,t.size):
    x_E[0][i] = U_mat[1][i-1]*np.cos(angle[2][i-1])*np.cos(angle[3][i-1])
    x_E[1][i] = U_mat[1][i-1]*np.cos(angle[2][i-1])*np.sin(angle[3][i-1])
    x_E[2][i] = U_mat[1][i-1]*np.sin(angle[2][i-1])
    for j in range(3,6):
        x_E[j][i] = x_E[j][i-1] + dif*x_E[j-3][i-1]
    D = rho*U_mat[1][i-1]**2*S/2*(C_d0+K*alpha**2)
    U_mat[0][i] = (-D+T*np.cos(alpha))/m-g*np.sin(angle[2][i-1])
    U_mat[1][i] = U_mat[1][i-1] + dif*U_mat[0][i-1]
    L = rho*U_mat[1][i-1]**2*S/2*C_lalpha*alpha
    angle[0][i] = ((L+T*np.sin(alpha))/m*np.cos(Phi[i])-g*np.cos(angle[2][i-1]))/U_mat[1][i-1]
    angle[1][i] = (L+T*np.sin(alpha))*np.sin(Phi[i])/m/np.cos(angle[2][i-1])/U_mat[1][i-1]
    angle[2][i] = angle[2][i-1] + dif*angle[0][i-1]
    angle[3][i] = angle[3][i-1] + dif*angle[1][i-1]

plt.figure()
plt.plot(t,U_mat[1])
plt.title("U")
plt.xlabel("t[s]")
plt.ylabel("U[m/s]")
plt.show()

plt.figure()
plt.plot(x_E[3],x_E[4])
plt.title("x-y")
plt.xlabel("x_e[m]")
plt.ylabel("y_e[m]")
plt.show()

plt.figure()
plt.plot(t,x_E[5])
plt.title("altitude")
plt.xlabel("t[s]")
plt.ylabel("h[m]")
plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x_E[3],x_E[4],x_E[5])
plt.title("x-y-z")
ax.set_xlabel('x_e[m]')
ax.set_ylabel('y_e[m]')
ax.set_zlabel('z_e[m]')
plt.show()