import numpy as np
import scipy.special as ss
import mpmath as mp

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
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N,1))
    a = np.zeros((N, M), np.complex128)


    if (scale == 1):
        for i in range(0, N):
            a[i,:] =  np.sqrt(np.pi*x/2) * ss.jve(nu[i]+0.5, x) 
    elif (scale == 0):
        for i in range(0, N):
            a[i,:] =  np.sqrt(np.pi*x/2) * ss.jv(nu[i]+0.5, x)
    else:
        print('incorrect input to ric_besselj')
    
    '''
    for i in range(0, len(nu)):
    	for j in range(0, len(x)):
            y = np.sqrt(np.pi * x[j] / 2) * mp.besselj(nu[i]+0.5, x[j]) * np.e**(-1*abs(np.imag(x[j])))
            a[i,j] = complex(y.real, y.imag)
    '''
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
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (len(nu), 1))
    
    #debugging print statements
    #print(x)
    #print(np.shape(x))
    #print(len(x))

    temp = np.matmul(np.ones((len(nu), 1)), x)

    if (flag == 1):
        #J = 0.5*( ric_besselj(nu-1, x) + (1/temp)*ric_besselj(nu,x) - ric_besselj(nu+1, x) )
        J0 = ric_besselj(nu-1, x)
        J1 = (1/temp)*ric_besselj(nu,x)
        J2 = ric_besselj(nu+1, x)
        J = 0.5*(J0+ J1 - J2)
    elif (flag == 2):
        J = 0.5 * ( ric_besselj_derivative(nu-1,x) + (1/temp)*ric_besselj_derivative(nu, x) \
                    - (1/(temp**2)) * ric_besselj(nu,x)  - ric_besselj_derivative(nu+1, x) )
    else:
        print('error: check third argument passed to ric_besselj_derivative (flag)')
    
    #removing all the zeros from inside the matrix...
    #f x*nu was 0, then it should become 0 after calculation
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
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N,1))
    a = np.zeros((N, M), np.complex128)

    if (scale == 1):
        for i in range(0, N):
            #print(nu, type(nu), type(nu[0]))
            #print(x, type(x), type(x[0]))
            a[i,:] =  np.sqrt(np.pi*x/2) * ss.yve(nu[i]+0.5, x) 
    elif (scale == 0):
        for i in range(0, N):
            a[i,:] =  np.sqrt(np.pi*x/2) * ss.yv(nu[i]+0.5, x)
    else:
        print('incorrect input to ric_besselj')
    
    # handling the case where x is zero because bessely is poorly defined 
    temp = np.matmul(np.ones((N,1)), x)
    a[temp == 0] = float('-inf')

    temp1 = np.matmul(nu, np.ones((1,M)))
    a[temp+temp1 == 0] = -1

    return a

def ric_bessely_derivative(nu, x, flag = 1):
    '''
        Y = ric_bessely_derivative(nu, x, flag) using the recursive
        relationship to calculate the first derivative of the
        Riccati-Neumann's function.

        The Riccati-Neumann's function is defined as
            Y_{nu}(x) = \sqrt{\frac{\pi x}{2}} Y_{nu+1/2}(x)

        Inputs:
        nu         order of the riccati Bessel's function. Must be a column vector.   
        x          row vector of size parameters at each frequency
        flag       1, first order derivative ; 2, second order derivative
    '''
    M = max(np.shape(x))
    N = max(np.shape(nu))
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N, 1))
    
    temp = np.matmul(np.ones((N, 1)), x)
    if (flag == 1):
        Y = 0.5 * (ric_bessely(nu-1, x) + (1.0/temp)* ric_bessely(nu,x) - ric_bessely(nu+1, x) )
        #Y0 = ric_bessely(nu-1, x);
        #Y1 = (1/temp)*ric_bessely(nu,   x);
        #Y2 = ric_bessely(nu+1, x);
        #Y = 0.5*(Y0+ Y1- Y2);
        #print("Y0:\n", Y0,"\nY1:\n", Y1, "\nY2:\n", Y2, "\nY:\n", Y)
    elif (flag == 2):
        Y = 0.5 * ( ric_bessely_derivative(nu-1,x) + (1/temp)*ric_bessely_derivative(nu, x) \
                    - (1/(temp**2)) * ric_bessely(nu,x)  - ric_bessely_derivative(nu+1, x) )
    else:
        print('error: third argument passed to ric_bessely_derivative must be 1 or 2')
    
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
            x   row vector of size parameters at each frequency
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

def spherical_besselj(nu,x, derivative = 0):
    '''
        Compute spherical bessel function of the first kind, using formula 
        given in: Abramowitz and Stegun, p. 437, 10.1.1
        Derivative calculated using recurrence relation in 
        (https://dlmf.nist.gov/10.51)

        Inputs:
        nu         order of the riccati Bessel's function. Must be a column vector.   
        x          row vector of size parameters at each frequency
        derivative       1, first order derivative ; 0, none.
    '''
    M = max(np.shape(x))
    N = max(np.shape(nu))
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N,1))
 
    j = ss.spherical_jn(nu, x, derivative)

    return j

def spherical_besselh(nu,x, derivative = 0):
    '''
        Compute spherical hankel function of the first kind, using formula 
        given in: Abramowitz and Stegun, p. 437, 10.1.1
        Derivative calculated using recurrence relation in 
        (https://dlmf.nist.gov/10.51)

        Inputs:
        nu         order of the riccati Bessel's function. Must be a column vector.   
        x          row vector of size parameters at each frequency
        derivative       1, first order derivative ; 0, none.
    '''
    M = max(np.shape(x))
    N = max(np.shape(nu))
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N,1))
 
    J = ss.spherical_jn(nu, x, derivative)
    Y = ss.spherical_yn(nu, x, derivative)
    H = J + (1.0j)*Y
    return H

if __name__ == "__main__":
    from DielectricMaterial import *
    radius = 0.5;
    sphere = DielectricMaterial(2.56,0.01)
    background = DielectricMaterial(1,0)
    sensor_location = [0,0,100];
    frequency = [1e6, 1e7, 1e8]

    x = []      #size parameter
    for f in frequency:
        x += [radius * sphere.getWaveNumber(f)]
    nu = np.arange(1,5,1)
    #using kzhu ricatti bessel
    result1 = ric_besselj_derivative(nu,x,1)

    M = max(np.shape(x))
    N = max(np.shape(nu))
    x = np.reshape(np.array(x), (1, M))
    nu = np.reshape(np.array(nu), (N,1))
    #using recurrence relations and scipy special
    result2 =   x*spherical_besselj(nu-1,x) - nu*spherical_besselj(nu,x)

    print(result1)
    print("")
    print(result2)

