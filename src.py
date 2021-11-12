import numpy as np
import scipy.special as ss

def cartToSph(x,y,z):
    '''
        [r, theta, phi] = cartToSph(x, y, z) converts the cartesian
        coordinate system to the spherical coordinate system according to
        the following definition:
        r       distance from the origin to the point in the interval 
                [0, \infty)
        theta   elevation angle measured between the positive z-axis and
                the vector in the interval [0, pi]
        phi     azimuth angle measured between the positive x-axis and
                the vector in the interval [0, 2*pi)
    '''

    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(np.sqrt(x**2 + y**2), z)
    phi = np.arctan2(y , x)

    return [r, theta, phi]

def sphToCart(r,theta,phi):
    '''
        [x,y,z] = sphToCart(r,theta,phi)
        for converting from spherical to cartesian coordinates
        r is radius, 
        theta is angle of elevation (0 at positive z axis, pi at negative z axis)
        phi is angle of azimuth (0 at positive x axis, increase counterclockwise)
    '''
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    return [x,y,z]

def get_legendre(n,m,x):
    '''
        Returns an array dimensions len(N) by 1 with the 
        value of the m-th degree term of the n-th order 
        associated legendre polynomial evaluated at x. 

        Inputs: 
            n: a sequence of integers  
            m: a single integer, for now.
            x: the argument to legenre polynomial
        Output: 
            P
    '''
    P = []
    for i in range(0,len(n)):
        # use scipy.special to computePmn(x)
        #into an m+1 x n+1 array for every value
        #0...m, 0...n.
        a,b = ss.lpmn(m,n[i],x)
        #select the value at the m,n of interest
        P.append(a[m,n[i]])
    return P

def norm(sensor_location):
    '''
        return the pythagorean distance from the sensor location
        to the origin. 
    '''
    x = sensor_location[0]
    y = sensor_location[1]
    z = sensor_location[2]
    return np.sqrt(x**2 + y**2 + z**2)



if __name__ == "__main__":
    sensor_location = [0,0,-2000]
    [r,theta,phi] = cartToSph(sensor_location)
    print(r,theta,phi)
	