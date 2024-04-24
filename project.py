import numpy as np
import matplotlib.pyplot as plt

#frequency of transmitter
fs = 1*10**9
#wave length 
wl = (3*10**8)/fs

#longatude latitude variables
long = 0
lat  = 0

#set location in ECEF
G = np.array([[6257.495860574267],
 [ 1103.3653518907781],
 [552.1839593406222]])

Q = 1**2*np.array([[2, 1, 1],
             [1, 2, 1],
             [1, 1, 2]])

Qinv = 1**2 * np.array([[ 0.75, -0.25, -0.25],
                [-0.25,  0.75, -0.25],
                [-0.25, -0.25, 0.75]])

#data matrix
theta = np.empty(shape = (41, 41)) 

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
    return x, y, z

#transmitter doppler shift
def Transmitter(G):
    Gf0 = 1/wl*((G - s1).T @ s1d)/(np.linalg.norm(G-s1))
    g1 = 1/wl*((G - s2).T @ s2d)/(np.linalg.norm(G-s2)) - Gf0
    g2 = 1/wl*((G - s3).T @ s3d)/(np.linalg.norm(G-s3)) - Gf0
    g3 = 1/wl*((G - s4).T @ s4d)/(np.linalg.norm(G-s4)) - Gf0
    return np.array([g1[0], g2[0], g3[0]])

#grid search doppler shift
def dopler_shift(u):
    f0 = 1/wl*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
    f1 = 1/wl*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - f0
    f2 = 1/wl*((u - s3).T @ s3d)/(np.linalg.norm(u-s3)) - f0
    f3 = 1/wl*((u - s4).T @ s4d)/(np.linalg.norm(u-s4)) - f0 
    return np.array([f1[0], f2[0], f3[0]])

#grid search calc
g = Transmitter(G)
for i in range(0, 41):
    for j in range(0, 41):
        u = np.array(ECEF_to_LLA(lat, long)).T
        f = dopler_shift(u)
        theta[j][i] = ((f - g).T @ Qinv @ (f - g))
        if i <= 5 and j<=10:
            print(lat, long)
            print((f-g).T@(f-g))
            print(theta[j][i])
        long += 1
    long = 0
    lat += 1

#print(np.min(theta))

#plotting stuff
feature_x = np.arange(0, 41, 1) 
feature_y = np.arange(0, 41, 1) 

[X, Y] = np.meshgrid(feature_x, feature_y)

contour = plt.contourf(X, Y, theta) 

plt.colorbar(contour, label='Z-values')

plt.show()