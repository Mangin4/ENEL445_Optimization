import numpy as np
import matplotlib.pyplot as plt

fs = 1*10**9      #frequency of transmitter
wl = (3*10**5)/fs #wave length 

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

def ECEF_to_LLA(x, y, z):
    a = 6378.1370 #km
    b = 6356.752314245 #km

    f = (a - b) / a
    f_inv = 1.0 / f
    e_sq = f * (2 - f)                       
    eps = e_sq / (1.0 - e_sq)
    p = np.sqrt(x * x + y * y)
    q = np.arctan2((z * a), (p * b))
    sin_q = np.sin(q)
    cos_q = np.cos(q)
    sin_q_3 = sin_q * sin_q * sin_q
    cos_q_3 = cos_q * cos_q * cos_q
    phi = np.arctan2((z + eps * b * sin_q_3), (p - e_sq * a * cos_q_3))
    lam = np.arctan2(y, x)
    v = a / np.sqrt(1.0 - e_sq * np.sin(phi) * np.sin(phi))
    h   = (p / np.cos(phi)) - v
    lat = np.degrees(phi)
    lon = np.degrees(lam)
    result_vec = np.array([[lat, lon, h]]).T
    return result_vec

#grid search doppler shift
def dopler_shift(lat, long, w):
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

def get_jac(x):
    u = LLA_to_ECEF((x[0])[0], (x[1])[0])
    f21d = 1/wl*(((s2d)/(np.linalg.norm(u-s2)))-(((u-s2).T @ s2d)*(u-s2))/((np.linalg.norm(u-s2))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    
    f31d = 1/wl*(((s3d)/(np.linalg.norm(u-s3)))-(((u-s3).T @ s3d)*(u-s3))/((np.linalg.norm(u-s3))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    
    f41d = 1/wl*(((s4d)/(np.linalg.norm(u-s4)))-(((u-s4).T @ s4d)*(u-s4))/((np.linalg.norm(u-s4))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    
    a = f21d[:, 0]
    b = f31d[:, 0]
    c = f41d[:, 0]
    J = np.array([a, b, c])
    lat = (x[0])[0]
    lon  = (x[1])[0]
    y = dopler_shift(5, 10, w=0)
    f = dopler_shift(lat, lon, w=0)
    r = y-f+J@x
    return J, r

def lev():
    x0 = np.array([[9],[10],[0]])
    mu0 = 1000
    ro = 10
    k = 0
    x = x0
    mu = mu0
    delta = 1
    j, r = get_jac(x)
    e = (np.linalg.norm(r))**2
    tol = 0.01
    while abs(delta) > tol:
        s = np.linalg.inv(j.T@j+mu*np.diag(np.diag(j.T@j)))@j.T@r
        js, rs = get_jac(x+s)
        es = (np.linalg.norm(rs))**2
        delta = es-e
        if delta < 0:
            x = x + s
            r, j, e = rs, js, es
            mu = mu/ro
        else:
            mu = mu*ro
        k+=1
        print(x)
    print(x, k)

def main():
    lev()

main()