import numpy as np

m = 7500
S = 27.9
g = 9.8
rho = 0.736423
U = 180
uhen = m*g*2/(rho*(U**2)*S)
print(uhen)

C_lalpha = 4.30
C_d0 = 0.0548
K = 3.02

dif = 1e-7

alpha = 0
while True:
    alpha += dif
    sahen = C_lalpha*alpha+C_d0*np.tan(alpha)+K*alpha**2*np.tan(alpha)
    if(sahen > uhen):
        break

print(alpha)
print(alpha/np.pi*180)
print(sahen)

T = (m*g-rho*U**2*S*C_lalpha*alpha/2)/np.sin(alpha)
print(T)

T= rho*U**2*S/2*(C_d0+K*alpha**2)/np.cos(alpha)
print(T)