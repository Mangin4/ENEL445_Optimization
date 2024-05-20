import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from numpy import ma
import time

MAX_LAT = 40
MAX_LONG = 40

lat_axis = np.linspace(0, MAX_LAT, 41)
long_axis = np.linspace(0, MAX_LONG, 41)

R0 = 6378.137 # Earth's radius
e = 0.081819198425 # Earth's eccentricity
fs = 1*10**9 # Source frequency
wl = (3*10**8)/fs # Source wavelength

p_source = np.array([[5], [10], [0]]) # Source location

Q = np.array([[2, 1, 1],
             [1, 2, 1],
             [1, 1, 2]])

Qinv = np.array([[ 0.75, -0.25, -0.25],
                [-0.25,  0.75, -0.25],
                [-0.25, -0.25, 0.75]])

#satalites
s1 = np.array([[7378.1, 0, 0]]).T
s1d = np.array([[0.0001, 4.4995, 5.3623]]).T
s2 = np.array([[7377.5, 100, 0]]).T
s2d = np.array([[-0.0671, 4.9493, 4.9497]]).T
s3 = np.array([[7377.5, -100, 0]]).T
s3d = np.array([[0.0610, 4.4991, 5.3623]]).T
s4 = np.array([[7377.5, 0, 100]]).T
s4d = np.array([[-0.0777, 4.0150, 5.7335]]).T

def convert(p):
    B = p[0, 0]
    L = p[1, 0]
    H = p[2, 0]

    Re = R0/np.sqrt((1-e**2*np.sin(np.pi*B/180)**2))
    
    x = (Re + H) * np.cos(np.pi*B/180) * np.cos(np.pi*L/180)
    y = (Re + H) * np.cos(np.pi*B/180) * np.sin(np.pi*L/180)
    z = ((1-e**2) * Re + H) * np.sin(np.pi*B/180)

    u = np.array([[x], [y], [z]])
    return u

u_source = convert(p_source)
f_ref = 1/wl * np.dot((u_source-s1).T, s1d)[0, 0] / np.linalg.norm(u_source - s1)

def getVec(p):
    u = convert(p)
    v_21 = 1/wl * np.dot((u-s2).T, s2d)[0, 0] / np.linalg.norm(u - s2) - f_ref
    v_31 = 1/wl * np.dot((u-s3).T, s3d)[0, 0] / np.linalg.norm(u - s3) - f_ref
    v_41 = 1/wl * np.dot((u-s4).T, s4d)[0, 0] / np.linalg.norm(u - s4) - f_ref
    result_vec = np.array([[v_21], [v_31], [v_41]])
    return result_vec

dummy = np.array([[5], [10], [0]])
#print(getVec(dummy))


def gridSearch():
    lat = 0
    long  = 0

    theta = np.empty(shape = (41, 41)) 
    f = getVec(p_source)
    for i in long_axis:
        for j in lat_axis:
            p_next = np.array([[lat], [long], [0]])
            g = getVec(p_next)
            theta[int(j)][int(i)] = np.log((f-g).T @ Qinv @ (f-g))[0, 0] # TODO THIS IS LOGGED BUT ITS WEIRD AT < 1 BECAUSE WILL BE NEGATIVE
            lat += 0.5
        lat = 0
        long += 0.5
            # WAS USING THIS TO TEST POINT 35, 35 long lat
            # if (j == 35 and i ==35):
            #     print(f"\nf: {f}")
            #     print(f"g: {g}")
            #     print(f"f-g: {f-g}")
            #     print(f".T {(f-g).T}")
            #     print(f"Qinv: {(f-g).T@Qinv}")
            #     print(f"WHOLE: {(f-g).T@Qinv@(f-g)}")
            #time.sleep(5)
    return theta

def main():

    theta = gridSearch()

    [X, Y] = np.meshgrid(lat_axis, long_axis)

    # fig, ax = plt.subplots()
    # cs = ax.contourf(X, Y, theta, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
    # cbar = fig.colorbar(cs)

    contour = plt.contourf(X, Y, theta) 

    plt.colorbar(contour, label='Z-values')

    # MORE TEST
    # print(np.max(theta))
    # print(np.min(theta))
    # print(theta[35][35])
    plt.show()
main()


