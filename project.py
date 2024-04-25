import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

#frequency of transmitter
fs = 1*10**9
#wave length 
wl = (3*10**8)/fs

#longatude latitude variables
long = 0
lat  = 0

Q = np.array([[2, 1, 1],
            [1, 2, 1],
            [1, 1, 2]])

Qinv = np.array([[ 0.75, -0.25, -0.25],
            [-0.25,  0.75, -0.25],
            [-0.25, -0.25, 0.75]])

#data matrix
theta = np.empty(shape = (40, 40)) 

#satalites
s1 = np.array([[7378.1, 0, 0]]).T
s1d = np.array([[0.0001, 4.4995, 5.3623]]).T
s2 = np.array([[7377.5, 100, 0]]).T
s2d = np.array([[-0.0671, 4.9493, 4.9497]]).T
s3 = np.array([[7377.5, -100, 0]]).T
s3d = np.array([[0.0610, 4.4991, 5.3623]]).T
s4 = np.array([[7377.5, 0, 100]]).T
s4d = np.array([[-0.0777, 4.0150, 5.7335]]).T

#Should convert from ECEF to LLA
def ECEF_to_LLA(lat, long):
    R0 = 6378137
    e = 0.081819198425
    lat = np.radians(lat)
    long = np.radians(long)
    Re = R0/(np.sqrt(1-(e)**2*(np.sin(lat))**2))
    x = (Re*np.cos(lat)*np.cos(long))
    y = (Re*np.cos(lat)*np.sin(long))
    z = ((1-e**2)*Re*np.sin(lat))
    u = np.array([[x], [y], [z]])
    return u

#grid search doppler shift
def dopler_shift(lat, long):
    u = ECEF_to_LLA(lat, long)
    f0 = 1/wl*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
    f1 = 1/wl*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - f0
    f2 = 1/wl*((u - s3).T @ s3d)/(np.linalg.norm(u-s3)) - f0
    f3 = 1/wl*((u - s4).T @ s4d)/(np.linalg.norm(u-s4)) - f0 
    f = np.array([f1[0], f2[0], f3[0]])
    return f

#grid search calc
f = dopler_shift(5, 10)
for i in range(0, 40):
    for j in range(0, 40):
        g = dopler_shift(lat, long)
        theta[j][i] = (f-g).T @ Qinv @ (f-g)
        lat += 1
    lat = 0
    long += 1

#plotting stuff
feature_x = np.arange(0, 40, 1) 
feature_y = np.arange(0, 40, 1) 

[X, Y] = np.meshgrid(feature_x, feature_y)

contour = plt.contourf(X, Y, theta, locator=ticker.LogLocator()) 

plt.colorbar(contour, label='Z-values')

plt.show()
