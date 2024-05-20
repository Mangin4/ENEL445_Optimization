import numpy as np
import matplotlib.pyplot as plt
import random
import time

fs = 1*10**9      #frequency of transmitter
wl = (3*10**5)/fs #wave length 
nVals = 4001
step = 0.01

#source location
u0 = np.array([[5], [10], [0]])

#Data matrix
theta = np.empty(shape = (nVals, nVals))

#Weighting function
Q = np.array([[2, 1, 1],
            [1, 2, 1],
            [1, 1, 2]])

Qinv = np.linalg.inv(Q)

#making the noise vector
n1 = random.gauss(0, 16)
n2 = random.gauss(0, 16)
n3 = random.gauss(0, 16)
n4 = random.gauss(0, 16)
nf = np.array([[n2-n1], [n3-n1], [n4-n1]])

#satalites
s1 = np.array([[7378.1, 0, 0]]).T
s2 = np.array([[7377.5, 100, 0]]).T
s3 = np.array([[7377.5, -100, 0]]).T
s4 = np.array([[7377.5, 0, 100]]).T

s1d = np.array([[0.0001, 4.4995, 5.3623]]).T
s2d = np.array([[-0.0671, 4.9493, 4.9497]]).T
s3d = np.array([[0.0610, 4.4991, 5.3623]]).T
s4d = np.array([[-0.0777, 4.0150, 5.7335]]).T

#Should convert from ECEF to LLA
def LLA_to_ECEF(lat, long):
    R0 = 6378.137
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
    u = LLA_to_ECEF(lat, long)
    f0 = 1/wl*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
    f1 = 1/wl*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - f0
    f2 = 1/wl*((u - s3).T @ s3d)/(np.linalg.norm(u-s3)) - f0
    f3 = 1/wl*((u - s4).T @ s4d)/(np.linalg.norm(u-s4)) - f0 
    f = np.array([f1[0], f2[0], f3[0]])
    return f

# def get_f(lat, long):
#     u = LLA_to_ECEF(lat, long)
#     f0 = ((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
#     f1 = ((u - s2).T @ s2d)/(np.linalg.norm(u-s2))
#     f2 = ((u - s3).T @ s3d)/(np.linalg.norm(u-s3))
#     f3 = ((u - s4).T @ s4d)/(np.linalg.norm(u-s4))
#     fs = (f0+f1+f2+f3)/4
#     return fs

#grid search calc
def map_gen():
    long = 0
    lat  = 0

    f = dopler_shift(5, 10)
    #get_f(5, 10)
    for i in range(0, nVals):
        for j in range(0, nVals):
            g = dopler_shift(lat, long)
            g += nf
            theta[j][i] = np.log((f - g).T @ Qinv @ (f - g))
            lat += step
        lat = 0
        long += step

def grid_search():
    coord = 0
    min = float('inf')
    for i in range(0,nVals):
         for j in range(0, nVals):
            if theta[j][i] < min: #grid search stuff
                min = theta[j][i]
                coord = np.array([[j], [i], [0]])*0.1
    print( coord, min)

#monte carlo simulation to find the accuracy of the data
def monte_carlo():
    sumation = 0
    L = 500
    for l in range(0, L):
        coord, min = grid_search()
        sumation += np.linalg.norm(coord - u0)**2
        
    RMSE = np.sqrt((1/L)*sumation)
    print(RMSE)

#plotting stuff
def graph():
    #coord, min = grid_search()
    tmin = np.min(theta)
    tmax = np.max(theta)
    levels = np.linspace(tmin, tmax, 100)
    feature_x = np.arange(0, 40.01, 0.01) 
    feature_y = np.arange(0, 40.01, 0.01) 

    [X, Y] = np.meshgrid(feature_x, feature_y)

    contour = plt.contourf(X, Y, theta, levels = levels) 

    plt.colorbar(contour, label='log(theta-values)')

    plt.xlabel("Longitude (Degrees)")
    plt.ylabel("Latitude (Degrees)")
    plt.title("FDOA plot")
    #print("Minmum is found at: ", coord.T)
    #plt.annotate('Min',xy=(10,5),xytext=(5,10),arrowprops={})
    
    plt.show()

def new_x(a, xc, xn):
    return xc + a*(xc-xn)
    
def main():
    start = time.time()
    map_gen()
    #grid_search()
    #monte_carlo()
    graph()
    #triang_thing()
    u = LLA_to_ECEF(5, 10)
    print(1/wl*(((s2d)/(np.linalg.norm(u-s2))**(1/2))-(((u-s2).T @ s2d)*(u-s2))/(2*(np.linalg.norm(u-s2))**3))
          -(1/wl*(((s1d)/(np.linalg.norm(u-s1))**(1/2))-(((u-s1).T @ s1d)*(u-s1))/(2*(np.linalg.norm(u-s1))**3))))
    end = time.time()
    timer = end - start
    print(timer/60)
main()  