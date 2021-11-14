# same functions as bessel.py 
# but using arbitrary precision library mpmath
# (https://mpmath.org/doc/0.19/index.html)

# Libraries that calculate bessel functions using fixed-precision 
# floating point variables may fail to evaluate besselj and bessely
# functions accuratelly. 
#
# For example, scipy.special 's jv(nu=0.5, x=9.4e14 + 1e-20j)
# returns (nan+nanj). MATLAB builtin besselj(0.5, 9.4e14 + 1e-20j)
# returns -1.7286e-08 + 1.6554e-24i. 
# For comparison, mpmath 's besselj(0.5, 9.4e14 + 1e-20j) returns
# mpc(real='-1.7286313...e-8', imag='-1.945349...e-28').
#  
# Clearly, using mpmath is preffered over scipy.special since the former 
# can compute more extreme inputs, which arise at high frequencies
# and high sphere conductivities. 
#
# To reap the benefits of arbitrary precision arithmetic, it 
# is sufficient to evaluate only the bessel functions using mpmath
# thereby avoiding scenarious where the result is (nan+nanj).

import numpy as np
import scipy.special as ss
import mpmath as mp

# set precision (default = 15 decimal places)
# for example, let's set to 50 decimal place precision
# this will be applied only to evaluate besselj,bessely functions
mp.dps = 50

def ric_besselj(nu,x,scale = 1):
    '''
        Implementation of riccatti bessel function 
        Inputs: 
            nu: column vector of integers from 1 to N inclusive
            x: rpw vector of M complex-valued arguments, M = len(frequency)
            scale: equals 1 by default, scales the output of bessel function 
                    by e**(-abs(imag(x))). if zero, does not scale. (not recommended)
        Output: 
            sqrt(pi*x/2)* J(nu+1/2, x)
            returned as an N by M array for each N values of nu and each M values of x
        Notes:    
            scale factors will cancel out in the end, and minimize numerical issues due to
            arguments with very large complex arguments. 
    '''
    M = max(np.shape(x))
    N = max(np.shape(nu))
    a = np.zeros((N, M), np.complex128)

    if (scale != 0 and scale != 1):
        print("incorrect argument ric_besselj (scale)")
        return
    
    for i in range(0, len(nu)):
    	for j in range(0, len(x)):
            if (scale == 1):
                y = np.sqrt(np.pi * x[j] / 2) * mp.besselj(nu[i]+0.5, x[j])* np.e**(-1*abs(np.imag(x[j])))
            elif (scale == 0):
                y = np.sqrt(np.pi * (x[j]) / 2) * mp.besselj(nu[i]+0.5, x[j])
            a[i,j] = complex(y.real, y.imag)

    return a

def ric_besselj_derivative(nu, x, flag = 1):
    '''
        translation of KZHU ric_besselj_derivative(nu,x,flag)
        Inputs:
            nu: order of Riccati-Bessel Function, integer sequence from 1:N as an array
            x: arguments to Ric-Bessel Function, as a vector with M elements, where M is 
                number of frequencies
            flag: 1 for first order derivative, 2 for second order derivative. 
    '''
    M = max(np.shape(x))
    M = max(np.shape(x))
    
    
    #debugging print statements
    #print(x)
    #print(np.shape(x))
    #print(len(x))

    temp = np.matmul(np.ones((len(nu), 1)), np.reshape(np.array(x), (1, M)))

    if (flag == 1):
        J = 0.5*( ric_besselj(nu-1, x) + (1/temp)*ric_besselj(nu,x) - ric_besselj(nu+1, x) )
    elif (flag == 2):
        J = 0.5 * ( ric_besselj_derivative(nu-1,x) + (1/temp)*ric_besselj_derivative(nu, x) \
                    - (1/(temp**2)) * ric_besselj(nu,x)  - ric_besselj_derivative(nu+1, x) )
    else:
        print('error: check third argument passed to ric_besselj_derivative (flag)')
    
    #removing all the zeros from inside the matrix...
    #f x*nu was 0, then it should become 0 after calculation
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (len(nu), 1))

    temp1 = np.matmul(np.ones((len(nu), 1)), x)
    J[temp1 == 0] = 0

    temp2 = np.matmul(nu, np.ones((1,len(x))))
    if (flag == 1):
        J[temp1+temp2 == 0] = 1;

    return J

def ric_bessely(nu, x, scale = 1):
    '''
    '''
    M = max(np.shape(x))
    N = max(np.shape(nu))
    a = np.zeros((N, M), np.complex128)


    if (scale != 0 and scale != 1):
        print("incorrect argument ric_bessely (scale)")
        return
    
    for i in range(0, len(nu)):
    	for j in range(0, len(x)):
            if (scale == 1):
                y = np.sqrt(np.pi * x[j] / 2) * mp.bessely(nu[i]+0.5, x[j]) * np.e**(-1*abs(np.imag(x[j])))
            elif (scale == 0):
                y = np.sqrt(np.pi * (x[j]) / 2) * mp.bessely(nu[i]+0.5, x[j])
            a[i,j] = complex(y.real, y.imag)
    
    # handling the case where x is zero because bessely is poorly defined 
    temp = np.matmul(np.ones((N,1)), np.reshape(np.array(x), (1, M)))
    a[temp == 0] = float('-inf')

    return a

def ric_bessely_derivative(nu, x, flag = 1):
    '''
        Y = ric_bessely_derivative(nu, x, flag) using the recursive
        relationship to calculate the first derivative of the
        Riccati-Neumann's function.

        The Riccati-Neumann's function is defined as
            Y_{nu}(x) = \sqrt{\frac{\pi x}{2}} Y_{nu+1/2}(x)

        Inputs:
        nu: order of Riccati-Bessel Function, integer sequence from 1:N as an array
            x: arguments to Ric-Bessel Function, as a vector with M elements, where M is 
                number of frequencies 
                Note: if x == 0, then a ValueError is raised (because bessely(nu,0)= -inf)
                      This should not happen because size parameters are non-zero
            flag: 1 for first order derivative, 2 for second order derivative. 
    '''
    M = max(np.shape(x))
    N = max(np.shape(nu))
    
    temp = np.matmul(np.ones((N, 1)), np.reshape(np.array(x), (1, M)))
    if (flag == 1):
        Y = 0.5 * (ric_bessely(nu-1, x) + (1.0/temp)* ric_bessely(nu,x) - ric_bessely(nu+1, x) )
    elif (flag == 2):
        Y = 0.5 * ( ric_bessely_derivative(nu-1,x) + (1/temp)*ric_bessely_derivative(nu, x) \
                    - (1/(temp**2)) * ric_bessely(nu,x)  - ric_bessely_derivative(nu+1, x) )
    else:
        print('error: third argument passed to ric_bessely_derivative must be 1 or 2')
    
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N, 1))

    temp2 = np.matmul(np.ones((N,1)), x)
    Y[temp2 == 0] = float('-inf')
    temp1 = np.matmul(nu, np.ones((1,M)))
    if (flag == 1):
        Y[temp1+temp2 == 0] = 1
    elif (flag == 2):
        Y[temp1+temp2 == 0] = -1
    
    return Y
    
def ric_besselh(nu,x, K):
    '''
        H = ric_besselh(nu, K, x) implement the Hankel function,
        which is defined as
            H_{nu}(x) = \sqrt{\frac{\pi x}{2}} H_{nu+1/2}(x)
        Inputs:
            nu  order of the spherical Bessel's function. Must be a column vector.
            x   Must be a row vector.
            K   1 for Hankel's function of the first kind; 2 for Hankel's
                function of the second kind.
    '''
    if (K == 1):
        H = ric_besselj(nu,x) + 1j*ric_bessely(nu,x)
    elif(K == 2):
        H = ric_besselj(nu,x) - 1j*ric_bessely(nu,x)
    else:
        print('error: third argument passed to ric_besselh must be 1 or 2')

    return H

def ric_besselh_derivative(nu, x, K, flag = 1):
    '''
        H = ric_besselh_derivative(nu, K, x) use the recursive relationship
        to calculate the first derivative of the Hankel function.
            H_{nu}(x) = \sqrt{\frac{\pi x}{2}} H_{nu+1/2}(x)

        Inputs:
            nu      order of the riccati-Hankel's function. Must be a column vector
           
            K = 1   if it is Hankel's function of the first kind; K=2 if it is 
                    Hankel's function of the second kind.  
            x       Must be a row evector
            flag    1 for the first order derivative; 2 for the second order derivative
    '''
    if (K == 1):
        H = ric_besselj_derivative(nu,x,flag) + 1j*ric_bessely_derivative(nu,x,flag)
    elif(K == 2):
        H = ric_besselj_derivative(nu,x,flag) - 1j*ric_bessely_derivative(nu,x,flag)
    else:
        print('error: argument K passed to ric_besselh_derivative must be 1 or 2')

    return H

if __name__ == "__main__":
    #debugging test script
    '''
    from DielectricMaterial import *
    radius = 0.5;
    sphere = DielectricMaterial(2.56,0.01)
    background = DielectricMaterial(1,0)
    sensor_location = [0,0,100];
    frequency = [1e6, 1e7, 1e8]

    x = []      #size parameter
    for f in frequency:
        x += [radius * sphere.getWaveNumber(f)]
    print(x)
    print(type(x[0]))

    nu = np.arange(1,5,1)
    print(nu)
    a = ric_besselj_derivative(nu,x,flag=1)
    print(a)
    '''