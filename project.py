import numpy as np
import matplotlib.pyplot as plt
import random
import time

fs = 1*10**9      #frequency of transmitter
wl = (3*10**5)/fs #wave length 
nVals = 401
step = 0.1

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

def gradient(u):
    derivitive = 1/wl*(((s2d)/(np.linalg.norm(u-s2)))-(((u-s2).T @ s2d)*(u-s2))/((np.linalg.norm(u-s2))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    return derivitive

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

    print(lat,lon,h)
    result_vec = np.array([[lat, lon, h]]).T
    return result_vec

#grid search doppler shift
def dopler_shift(lat, long, w):
    u = LLA_to_ECEF(lat, long)
    f0 = 1/w*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))
    f1 = 1/w*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - f0
    f2 = 1/w*((u - s3).T @ s3d)/(np.linalg.norm(u-s3)) - f0
    f3 = 1/w*((u - s4).T @ s4d)/(np.linalg.norm(u-s4)) - f0 
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
    freq = get_f(5, 10)
    w = (3*10**5)/freq
    f = dopler_shift(5, 10, w)
    for i in range(0, nVals):
        for j in range(0, nVals):
            g = dopler_shift(lat, long, w)
            #g += nf
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
    print(coord, min)

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
    feature_x = np.arange(0, 40+step, step) 
    feature_y = np.arange(0, 40+step, step) 

    [X, Y] = np.meshgrid(feature_x, feature_y)

    contour = plt.contourf(X, Y, theta, levels = levels) 

    plt.colorbar(contour, label='log(theta-values)')

    plt.xlabel("Longitude (Degrees)")
    plt.ylabel("Latitude (Degrees)")
    plt.title("FDOA plot")
    #print("Minmum is found at: ", coord.T)
    #plt.annotate('Min',xy=(10,5),xytext=(5,10),arrowprops={})
    
    plt.show()
    

def finite_diff(starting_pos, step_size, direction, wavelength):
    u = starting_pos #LLA_to_ECEF(starting_pos[0], starting_pos[1])
    h = step_size
    e = direction.T
    w = wavelength

    deriv = ((1/w*((u + h*e - s2).T @ s2d)/(np.linalg.norm(u + h*e - s2)) - 1/w*((u + h*e - s1).T @ s1d)/(np.linalg.norm(u + h*e - s1))) 
    - (1/w*((u - s2).T @ s2d)/(np.linalg.norm(u-s2)) - 1/w*((u - s1).T @ s1d)/(np.linalg.norm(u-s1))))/h

    return deriv[0][0]


def find():
    tol = 0.1
    pgrad = [0,0,0]
    grad = 0
    k = 0
    ppk = 0
    pk = 0
    bk = 0
    ak = 500
    pak = 0
    x = np.array([[6067.8,1625.9,1100.2]]).T

    grad = gradient(x)
    #np.linalg.norm(grad, ord = np.inf) > tol
    
    while k< 10:
        grad = gradient(x)
        ECEF_to_LLA(x[0], x[1], x[2])
        if k % 2 or k == 0:
            pk = -((grad)/np.linalg.norm(grad))
        else:
            bk = (np.dot(grad.T,grad)/np.dot(pgrad.T,pgrad))
            pk = -(grad)/np.linalg.norm(grad) + bk*ppk
        if k > 0:
            ak = pak*(np.dot(pgrad.T,ppk)/np.dot(grad.T,pk))
        x = x + ak*pk
        pgrad = np.linalg.norm(grad)
        ppk = pk
        pak = ak
        print(f"ppk {(ppk)}")
        print(f"pk {(pk)}")
        print(f"ak {(ak)}")
        print(f"bk {(bk)}")
        k += 1
    
    print(ECEF_to_LLA(x[0], x[1], x[2]))

# Function to compute Jacobian at position u and return 3x3 (includes altitude, H)
def get_jac(u):
    u = LLA_to_ECEF((u[0])[0], (u[1])[0])
    f21d = 1/wl*(((s2d)/(np.linalg.norm(u-s2)))-(((u-s2).T @ s2d)*(u-s2))/((np.linalg.norm(u-s2))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    
    f31d = 1/wl*(((s3d)/(np.linalg.norm(u-s3)))-(((u-s3).T @ s3d)*(u-s3))/((np.linalg.norm(u-s3))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    
    f41d = 1/wl*(((s4d)/(np.linalg.norm(u-s4)))-(((u-s4).T @ s4d)*(u-s4))/((np.linalg.norm(u-s4))**3))
    -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3)))
    
    a = f21d[:, 0]
    b = f31d[:, 0]
    c = f41d[:, 0]
    jacobian = np.array([a, b, c])
    return jacobian

# Should be correct new step but is far too large. Eventually H will get large and turn matrix singular.
def new_step(curr_pos, step_size):
    # [Bnew, Lnew] = (J.T @ J)^(-1) @ (J.T @ g), g = y - f(Bcurr, Lcurr)+J(Bcurr, Lcurr)[Bcurr, Lcurr]
    # Need this order of operation to get 2x1 result
    mu = step_size
    u = curr_pos # Set current position - use this for function arguments (has all 3 components)
    #u_2x1 = np.array([curr_pos[0], curr_pos[1]]) # Use this for the pos vector that is multiplied.
    B = (curr_pos[0])[0] # Get Longitude
    L = (curr_pos[1])[0] # Get Lattitude
    y = dopler_shift(5, 10, wl) # Doppler shift of target (will use estimate here)
    f = dopler_shift(B, L, wl) # Doppler shift of guess
    J = get_jac(u) # Get jacobian at current position (uses 3x1 vector in get_jac func)
    g = (y-f)+J @ u # Compute g used in finding step equation 
    D = np.diag(np.diag(mu * (J.T @ J)))
    J_inv = np.linalg.inv(J.T @ J + D) # Get inverse part of equation for new step
    u_next = -1*(J_inv @ (J.T @ g)) # Compute new step
    return u_next


# Not sure if this is how to evaluate the objective function.
def evaluate(curr_pos):
    y = dopler_shift(5, 10, wl)
    f = dopler_shift((curr_pos[0])[0], (curr_pos[1])[0], wl)
    res = (y-f).T @ (y-f)
    return res

# Pitiful attempt at implementing shit
def levenberg(starting_pos):
    tol = 500
    step = 0.1
    k = 0
    curr_pos = starting_pos
    #res = (evaluate(curr_pos)[0])[0]
    while (k < 20): # Think this is fine?
        new_pos = new_step(curr_pos, step)
        #print(f"New pos {new_pos}")
        new_pos = ECEF_to_LLA((new_pos[0])[0], (new_pos[1])[0], (new_pos[2])[0])
        new_pos = curr_pos - new_pos
        #res = evaluate(new_pos)
        curr_pos = new_pos
        print(new_pos)
        k+=1
    print(ECEF_to_LLA((new_pos[0])[0], (new_pos[1])[0], (new_pos[2])[0]))
    print("Converged")

def grad_descent(starting_pos):
    curr_pos = LLA_to_ECEF((starting_pos[0])[0], (starting_pos[1])[0])

    for i in range(0, 20):
        jac = get_jac(curr_pos)
        #new_pos = ECEF_to_LLA((new_pos[0])[0], (new_pos[1])[0], (new_pos[2])[0])
        new_pos = curr_pos - 200 * (jac @ curr_pos)
        curr_pos = new_pos
        print(curr_pos)
    print("CONVERGED")
    print(ECEF_to_LLA((new_pos[0])[0], (new_pos[1])[0], (new_pos[2])[0]))

def main():
    #start = time.time()
    #map_gen()
    #find()
    #grid_search()
    #monte_carlo()
    #graph()
    #triang_thing()
    # u = LLA_to_ECEF(5, 10)
    # print(1/wl*(((s2d)/(np.linalg.norm(u-s2)))-(((u-s2).T @ s2d)*(u-s2))/((np.linalg.norm(u-s2))**3))
    #        -(1/wl*(((s1d)/(np.linalg.norm(u-s1)))-(((u-s1).T @ s1d)*(u-s1))/((np.linalg.norm(u-s1))**3))))
    #end = time.time()
    #timer = end - start
    #print(timer/60)

    # Gradient Approximation
    # B_deriv = finite_diff(u, 1, np.array([[1, 0, 0]]), wl)
    # L_deriv = finite_diff(u, 1, np.array([[0, 1, 0]]), wl)
    # H_deriv = finite_diff(u, 1, np.array([[0, 0, 1]]), wl)
 
    # Test if converges when super close to point (it doesn't). It goes unstable af.
    start_pos = np.array([[20, 20, 0]]).T
    #levenberg(start_pos)
    grad_descent(start_pos)
    # a = new_step(start_pos, 0.01)
    # print(a)
    # a = ECEF_to_LLA((a[0])[0], (a[1])[0], (a[2])[0])


main()  