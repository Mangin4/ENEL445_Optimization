import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal
import plotly.graph_objects as go
from scipy.stats import multivariate_normal

fs = 1*10**9
wl = (3*10**8)/fs

STD_DEV = [1, 2, 4, 8, 16]

G = np.array([[-1513.22940634],
 [ -981.11867003],
 [-6055.80033183]])

Q = 1**2*np.array([[2, 1, 1],
             [1, 2, 1],
             [1, 1, 2]])



Qinv = 1**2 * np.array([[ 0.75, -0.25, -0.25],
                [-0.25,  0.75, -0.25],
                [-0.25, -0.25, 0.75]])

theta = np.empty(shape = (40, 40))

s1 = np.array([[7378.1, 0, 0]]).T
s1d = np.array([[0.0001, 4.4995, 5.3623]]).T
s2 = np.array([[7377.5, 100, 0]]).T
s2d = np.array([[-0.0671, 4.9493, 4.9497]]).T
s3 = np.array([[7377.5, -100, 0]]).T
s3d = np.array([[0.0610, 4.4991, 5.3623]]).T
s4 = np.array([[7377.5, 0, 100]]).T
s4d = np.array([[-0.0777, 4.0150, 5.7335]]).T

#Should convert from ECEF to LLA
def ECEF_to_LLA(p):
    R0 = 6378.137
    e = 0.081819198425
    Re = R0/(np.sqrt(1-e**2*np.sin(p[0, 0])))
    x = (Re*np.cos(p[0, 0])*np.cos(p[0, 1]))
    y = (Re*np.cos(p[0, 0])*np.sin(p[0, 1]))
    z = ((1-e**2)*Re*np.sin(p[0,0]))
    return x, y, z


def dopler_shift(u):
    f0 = 1/wl*(np.dot((u - s1).T, s1d))/(np.linalg.norm(u-s1))
    f1 = 1/wl*(np.dot((u - s2).T, s2d))/(np.linalg.norm(u-s2)) - f0 + np.random.rand(1)
    f2 = 1/wl*(np.dot((u - s3).T, s3d))/(np.linalg.norm(u-s3)) - f0 + np.random.rand(1)
    f3 = 1/wl*(np.dot((u - s4).T, s4d))/(np.linalg.norm(u-s4)) - f0 + np.random.rand(1)
    g1 = 1/wl*(np.dot((u - s2).T, s2d))/(np.linalg.norm(u-s2)) - f0
    g2 = 1/wl*(np.dot((u - s3).T, s3d))/(np.linalg.norm(u-s3)) - f0
    g3 = 1/wl*(np.dot((u - s4).T, s4d))/(np.linalg.norm(u-s4)) - f0
    return np.array([f1[0], f2[0], f3[0]]), np.array([g1[0], g2[0], g3[0]])

p = np.array([[0, 0, 0]])

for i in range(0, 39):
    for j in range(0, 39):
        u = np.array([ECEF_to_LLA(p)]).T
        f, g = dopler_shift(u)
        #f += np.random.rand(3,1)
        theta[j][i] = ((f - g).T @  Qinv) @ (f - g)
        p[0][0] += 1
    p[0][1] += 1

x = np.arange(0, 40, 1)
y = np.arange(0, 40, 1)

fig = go.Figure(data =
    go.Contour(z = theta,
    x = x,
    y = y))
fig.show()


