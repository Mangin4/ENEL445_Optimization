import numpy as np
import matplotlib.pyplot as plt
import random
import time

nVals = 401
step = 0.1

R0 = 6378.137
e = 0.081819198425

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
    lat = np.radians(lat)
    long = np.radians(long)
    Re = R0/(np.sqrt(1-(e)**2*(np.sin(lat))**2))
    x = (Re*np.cos(lat)*np.cos(long))
    y = (Re*np.cos(lat)*np.sin(long))
    z = ((1-e**2)*Re*np.sin(lat))
    u = np.array([[x], [y], [z]])
    return u

u0 = LLA_to_ECEF(5, 10)[:2,0]

#grid search doppler shift
def dopler_shift(lat, long, wl):
    u = LLA_to_ECEF(lat, long)
    f0 = 1/wl*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
    f1 = 1/wl*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - f0
    f2 = 1/wl*((u - s3).T @ s3d)/(np.linalg.norm(u-s3)) - f0
    f3 = 1/wl*((u - s4).T @ s4d)/(np.linalg.norm(u-s4)) - f0 
    f = np.array([f1[0], f2[0], f3[0]])
    return f

def get_f(lat, long):
    u = LLA_to_ECEF(lat, long)
    f0 = ((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
    f1 = ((u - s2).T @ s2d)/(np.linalg.norm(u-s2))
    f2 = ((u - s3).T @ s3d)/(np.linalg.norm(u-s3)) 
    f3 = ((u - s4).T @ s4d)/(np.linalg.norm(u-s4))
    fs = (f0+f1+f2+f3)/4
    return fs

#grid search calc
def map_gen():
    long = 0
    lat  = 0
    coord = 0
    min = float('inf')
    fe = get_f(5, 10)/(3*10**5)
    f = dopler_shift(5, 10, fe)
    for i in range(0, nVals):
        for j in range(0, nVals):
            g = dopler_shift(lat, long, fe)
            g += nf
            theta[j][i] = np.log((f - g).T @ Qinv @ (f - g))
            if theta[j][i] < min: #grid search stuff
                min = theta[j][i]
                coord = np.array([[j], [i], [0]])*0.1
            lat += step
        lat = 0
        long += step
    return coord, min

def finite_diff(starting_pos, step_size, direction, wavelength):
    u = starting_pos #LLA_to_ECEF(starting_pos[0], starting_pos[1])
    h = step_size
    e = direction.T
    w = wavelength
    deriv = ((1/w*((u + h*e - s2).T @ s2d)/(np.linalg.norm(u + h*e - s2)) - 1/w*((u + h*e - s1).T @ s1d)/(np.linalg.norm(u + h*e - s1))) 
    - (1/w*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - 1/w*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))))/h

    return deriv[0][0]


def complex_step(pos, step_size, direction, w):
    # w is estimated source wavelength
    h = step_size # Set to 10^(-200)
    e = direction.T
    u = pos + 1j*h*e # Perturb function with complex step increment
    u = np.array(u, dtype=complex)
    deriv = (np.imag(1/w*((u - s2).T @ s2d)/(np.linalg.norm(u - s2)) - 1/w*((u - s1).T @ s1d)/(np.linalg.norm(u - s1))))/h 
    return deriv[0][0]

def pos_jac(lat, lon):
    lat = np.radians(lat)
    lon = np.radians(lon)
    DRe = R0*(1-e**2*np.sin(lat)**2)**(-2/3)*e**2*np.sin(lat)*np.cos(lat)
    Re = R0/(np.sqrt(1-(e)**2*(np.sin(lat))**2))
    Bx = np.cos(lon)*(DRe*np.cos(lat)-Re*np.cos(lat))
    Lx = -Re*np.cos(lat)*np.sin(lon)
    By = np.sin(lon)*(DRe*np.cos(lat)-Re*np.sin(lat))
    Ly = Re*np.cos(lat)*np.cos(lon)
    Bz = (1-e**2)*(DRe*np.sin(lat)+Re*np.cos(lat))
    Lz = 0
    pos = np.array([[Bx, Lx],
                    [By, Ly],
                    [Bz, Lz]])
    return pos

def get_jac(x, wl):
    u = LLA_to_ECEF((x[0])[0], (x[1])[0])
    pos = pos_jac((x[0])[0], (x[1])[0])
    f21d = np.pi*(1/wl*(((s2d)/(np.linalg.norm(u-s2)))-(((u-s2).T @ s2d)*(u-s2))/((np.linalg.norm(u-s2))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3))))/180
    
    f31d = np.pi*(1/wl*(((s3d)/(np.linalg.norm(u-s3)))-(((u-s3).T @ s3d)*(u-s3))/((np.linalg.norm(u-s3))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3))))/180
    
    f41d = np.pi*(1/wl*(((s4d)/(np.linalg.norm(u-s4)))-(((u-s4).T @ s4d)*(u-s4))/((np.linalg.norm(u-s4))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3))))/180
    
    a = f21d[:, 0]
    b = f31d[:, 0]
    c = f41d[:, 0]
    J = np.array([a, b, c])@pos
    lat = (x[0])[0]
    lon  = (x[1])[0]
    y = dopler_shift(5, 10, wl)
    f = dopler_shift(lat, lon, wl)
    r = y-f+nf
    return J, r

def lev():
    x0 = np.array([[9],[10]])
    fe = get_f(5, 10)/(3*10**5)
    xs = np.empty(shape = (0))
    ys = np.empty(shape = (0))
    flag = 0
    mu0 = 0.1
    ro = 2
    k = 0
    x = x0
    mu = mu0
    delta = 1
    j, r = get_jac(x, fe)
    e = (np.linalg.norm(r))**2
    tol = 10**-9
    xs = np.append(xs, x[0])
    ys = np.append(ys, x[1])
    
    while abs(delta) > tol:
        s = np.linalg.inv(j.T@Qinv@j+mu*np.diag(np.diag(j.T@Qinv@j)))@j.T@Qinv@r
        if k < 2:
            s /=2
        js, rs = get_jac(x+s, fe)
        es = (np.linalg.norm(rs))**2
        delta = es-e
        if delta < 0:
            x = x + s
            xs = np.append(xs, x[0])
            ys = np.append(ys, x[1])
            r, j, e = rs, js, es
            if k % 5:
                mu = mu/ro
        else:
            mu = mu*ro
        if x[0] < 0 or x[0] > 40 or x[1] < 0 or x[1] > 40:
            flag = 1
            break
        k+=1
    if flag:
        print("not converged")
    # else:
    #     print(x, k)
    return x[:2,0], xs, ys

#plotting stuff
def graph(xs, ys):
    tmin = np.min(theta)
    tmax = np.max(theta)
    levels = np.linspace(tmin, tmax, 100)
    feature_x = np.arange(0, 40+step, step) 
    feature_y = np.arange(0, 40+step, step) 

    [X, Y] = np.meshgrid(feature_x, feature_y)

    contour = plt.contourf(X, Y, theta, levels = levels) 
    path = plt.plot(ys, xs) 

    plt.colorbar(contour, label='log(theta-values)')

    plt.xlabel("Longitude (Degrees)")
    plt.ylabel("Latitude (Degrees)")
    plt.title("FDOA plot")
    #print("Minmum is found at: ", coord.T)
    #plt.annotate('Min',xy=(10,5),xytext=(5,10),arrowprops={})
    
    plt.show()

#monte carlo simulation to find the accuracy of the data
def monte_carlo():
    sumation = 0
    L = 10**4
    for l in range(0, L):
        #val, min = map_gen()
        val, xs, ys = lev()
        print(val[1])
        val = LLA_to_ECEF(val[0], val[1])[:2,0]
        sumation += np.linalg.norm(val - u0)**2
        
    RMSE = np.sqrt((1/L)*sumation)
    print(RMSE)

def main():
    #map_gen()
    #val, xs, ys = lev()
    #monte_carlo()
    #graph(xs, ys)

    fe = get_f(5, 10)/(3*10**5)

    # Test vectors for derivative approximation and verification
    vec = np.array([[5], [10]])
    u = np.array([[5, 10, 0]]).T

    # Analytical df21/dx solution
    df21x = (1/fe*(((s2d)/(np.linalg.norm(u-s2)))-(((u-s2).T @ s2d)*(u-s2))/((np.linalg.norm(u-s2))**3))
    -(1/fe*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3))))
    
    print("Analytical Derivatives:")
    print(f"B: {df21x[0, 0]}")
    print(f"L: {df21x[1, 0]}")
    print(f"H: {df21x[2, 0]}")

    # Finite Difference Approximation
    B_deriv_fd = finite_diff(u, 10**(-8), np.array([[1, 0, 0]]), fe)
    L_deriv_fd = finite_diff(u, 10**(-8), np.array([[0, 1, 0]]), fe)
    H_deriv_fd = finite_diff(u, 10**(-8), np.array([[0, 0, 1]]), fe)

    print("Finite Difference Derivatives:")
    print(f"B: {B_deriv_fd}")
    print(f"L: {L_deriv_fd}")
    print(f"H: {H_deriv_fd}")

    # Complex Step Approximation
    B_deriv_cs = complex_step(u, 10**(-200), np.array([[1, 0, 0]]), fe)
    L_deriv_cs = complex_step(u, 10**(-200), np.array([[0, 1, 0]]), fe)
    H_deriv_cs = complex_step(u, 10**(-200), np.array([[0, 0, 1]]), fe)

    print(f"Complex Step Derivatives:")
    print(f"B: {B_deriv_cs}")
    print(f"L: {L_deriv_cs}")
    print(f"H: {H_deriv_cs}")

main()  