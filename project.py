import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, multivariate_normal

f0 = 1*10**9
wl = (3*10**8)/f0

STD_DEV = [1, 2, 4, 8, 16]

Q = np.array([[2, 1, 1],
             [1, 2, 1],
             [1, 1, 2]])

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
    f = 1/wl*(np.dot((u - s2).T, s2d))/(np.linalg.norm(u-s2))-1/wl*(np.dot((u - s1).T, s1d))/(np.linalg.norm(u-s1)) + np.random.rand(1)
    g = 1/wl*(np.dot((u - s2).T, s2d))/(np.linalg.norm(u-s2))-1/wl*(np.dot((u - s1).T, s1d))/(np.linalg.norm(u-s1))
    return f, g

p = np.array([[5,10, 0]])
u = np.array([ECEF_to_LLA(p)]).T
f, g = dopler_shift(u)
print(f)




