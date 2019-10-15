from numpy import *
import numpy.linalg as LA
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import *

a = 235*1e-3
b = 182*1e-3
h = 1.0*1e-3
E = 67*1e9
mu = 0.326
rho = 2.66*1e3
D = E*h**3/(12*(1-mu**2))
p = 5
q = 7
epsilon = np.empty(p)
epsilon[0] = 1.8751041
epsilon[1] = 4.6940911
epsilon[2] = 7.8547574
epsilon[3] = 10.9955407
epsilon[4] = 14.1371684

chi = np.empty(q)
chi[0] = 0
chi[1] = 0
chi[2] = 4.7300408
chi[3] = 7.8532046
chi[4] = 10.9956078
chi[5] = 14.1371655
chi[6] = 17.2787596

G = np.empty((p,q))
for i in range(p):
    for j in range(q):
        G[i][j] = b/a*epsilon[i]**4 + (a/b)**3*chi[j]**4
E = np.empty((p,p))
H = np.empty((p,p))
F = np.zeros((q,q))
K = np.zeros((q,q))

E[0][0] = 0.85824
E[0][1] = -11.74322
E[0][2] = 27.45315
E[0][3] = -37.39025
E[0][4] = 51.95662
E[1][0] = 1.87385
E[1][1] = -13.29425
E[1][2] = -9.04222
E[1][3] = 30.40219
E[1][4] = -33.70907
E[2][0] = 1.56451
E[2][1] = 3.22933
E[2][2] = -45.90423
E[2][3] = -8.33537
E[2][4] = 36.38656
E[3][0] = 1.08737
E[3][1] = 5.54065
E[3][2] = 4.25360
E[3][3] = -98.91821
E[3][4] = -7.82895
E[4][0] = 0.91404
E[4][1] = 5.71642
E[4][2] = 11.23264
E[4][3] = 4.73605
E[4][4] = -171.58466

H[0][0] = 4.64778
H[0][1] = -7.37987
H[0][2] = 3.94151
H[0][3] = -6.59339
H[0][4] = 4.59198
H[1][1] = 32.41735
H[1][2] = -22.35243
H[1][3] = 13.58425
H[1][4] = -22.83952
H[2][2] = 77.29889
H[2][3] = -35.64827
H[2][4] = 20.16203
H[3][3] = 142.90185
H[3][4] = -48.71964
for i in range(p):
    for j in range(i):
        H[i][j] = H[j][i]

F[0][2] = 18.58910
F[0][4] = 43.98096
F[0][6] = 69.11504
F[1][3] = 40.59448
F[1][5] = 84.08889
F[2][2] = -12.30262
F[2][4] = 52.58440
F[2][6] = 101.62255
F[3][3] = -46.05012
F[3][5] = 55.50868
F[4][2] = 1.80069
F[4][4] = -98.90480
F[4][6] = 60.12891
F[5][3] = 5.28566
F[5][5] = -171.58566
F[6][2] = 0.57069
F[6][4] = 9.86075
F[6][6] = -263.99798

K[1][1] = 12.00000
K[1][3] = 13.85641
K[1][5] = 13.85641
K[2][2] = 49.48082
K[2][4] = 35.37751
K[2][6] = 36.60752
K[3][3] = 108.92459
K[3][5] = 57.58881
K[4][4] = 186.86671
K[4][6] = 78.10116
K[5][5] = 284.68314
F[6][6] = 402.22805
for i in range(1,q):
    for j in range(1,i):
        K[i][j] = K[j][i]

C = np.empty((p,p,q,q))
for i in range(C.shape[0]):
    for m in range(C.shape[1]):
        for k in range(C.shape[2]):
            for n in range(C.shape[3]):
                C[i][m][k][n] = mu*a/b*(E[i][m]*F[n][k]+E[m][i]*F[k][n])+2*(1-mu)*a/b*H[i][m]*K[k][n]

J = np.zeros((p*q,p*q))
N = np.empty((p*q,p*q))
for i in range(p*q):
    J[i][i] = G[i//q][i%q]
for i in range(p*q):
    for j in range(p*q):
        N[i][j] = C[i//q][j//q] [i%q][j%q]
M = J+N
lmd_list, v = LA.eig(M)
v = v.T

lmd_vector_list = []
for i in range(len(lmd_list)):
    lmd_vector_list.append((lmd_list[i],v[i]))

lmd_vector_list.sort(key = lambda x:x[0])

omega_list = []
for i in range(len(lmd_vector_list)):
    omega = sqrt(lmd_vector_list[i][0]*D/(rho*h*a**3*b))
    omega_list.append(omega)

a *= 1e3
b *= 1e3

def Xm(x,m):
    alpha = (sinh(epsilon[m])-sin(epsilon[m]))/(cosh(epsilon[m])+cos(epsilon[m]))
    return cosh(epsilon[m]*x/a)-cos(epsilon[m]*x/a)-alpha*(sinh(epsilon[m]*x/a)-sin(epsilon[m]*x/a))

def Yn(y,n):
    if n == 0:
        return 1.0
    elif n == 1:
        return sqrt(3)*(1-2*y/b)
    else:
        beta = (sinh(chi[n])+sin(chi[n]))/(cosh(chi[n])-cos(chi[n]))
        return cosh(chi[n]*y/b)+cos(chi[n]*y/b)-beta*(sinh(chi[n]*y/b)+sin(chi[n]*y/b))

def W(x,y,i):
    W = 0
    for m in range(p):
        for n in range(q):
            W += lmd_vector_list[i][1][m*q+n]*Xm(x,m)*Yn(y,n)
    return W

x = arange(0,a,a/100)
y = arange(0,b,b/100)
X,Y = meshgrid(x,y)
for i in range(p*q):
    Z = W(X,Y,i)
    fig = figure()
    ax = fig.add_subplot(111, projection = "3d")
    surf = ax.plot_surface(X,Y,Z, cmap = "winter")
    ax.contourf(X, Y, Z, zdir='z', offset=-10)
    ax.contour(X, Y, Z, zdir='x', offset=-30)
    ax.contour(X, Y, Z, zdir='y', offset=b+30)
    ax.set_xlabel(r"$x$ [mm]")
    ax.set_ylabel(r"$y$ [mm]")
    ax.set_zlabel('$W$')
    ax.set_xlim3d(-30,a)
    ax.set_ylim3d(0, b+30)
    ax.set_zlim3d(-10, 10)
    ax.set_title(r"mode%d:  $\omega$=%s[Hz]"%(i+1,omega_list[i]/(2*pi)))
    colorbar(surf)
    subplots_adjust(left=-0.08, right=1.03, bottom=0.03, top=0.97)
    show(surf)
