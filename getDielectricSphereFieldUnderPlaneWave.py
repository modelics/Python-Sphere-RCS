import numpy as np
from DielectricMaterial import DielectricMaterial as DM
from src import *
from bessel import *
from getNMax import *

def getDielectricSphereFieldUnderPlaneWave(radius, sphere, background, sensor_location, frequency):
    '''
        Calculate the field as a plane wave is scattered by a dielectric
        sphere centered at the origin. The incident plane wave is
        polarized in the +x direction propagating in the +z
        direction.

        Inputs:
        radius             Scalar to denote the radius of the sphere (m)
        sphere             Object of DielectricMaterial
        background         Object of DielectricMaterial
        sensor_location    3x1 vector (m)
        frequency          Nx1 vector (Hz)

        Outputs:
        E_r, E_phi, E_theta     Nx1 vector (V/m)
        H_r, H_phi, H_theta     Nx1 vector (A/m)           
       
    '''
    

    EPS_0   = 8.8541878176*1e-12
    MU_0    = 4*np.pi*1e-7

    #all these values are of class 'numpy.ndarray'
    eta_m   = DM.getIntrinsicImpedance(background, frequency)
    k_m     = DM.getWaveNumber(background, frequency)
    mu_m    = DM.getComplexPermeability(background, frequency)*MU_0
    eps_m   = DM.getComplexPermittivity(background, frequency)*EPS_0
    
    eta_s   = DM.getIntrinsicImpedance(sphere, frequency)
    k_s     = DM.getWaveNumber(sphere, frequency)
    mu_s    = DM.getComplexPermeability(sphere, frequency)*MU_0
    eps_s   = DM.getComplexPermittivity(sphere, frequency)*EPS_0

    #debugging prints
    ''' 
    print(eta_m)
    print(k_m)
    print(mu_m)
    print(eps_m)

    print(eta_s)
    print(k_s)
    print(mu_s)
    print(eps_s)
    '''

    # number of mie terms to evaluate
    #based off wiscombe recommendation
    N = getNMax(radius, sphere, background, frequency)
    print("NMax: ", N)

    #need to convert frequency into numpy array
    if (type(frequency) == int or type(frequency) == float):
        frequency = [frequency]
    if (type(frequency) == list or type(frequency) == np.ndarray):
        frequency = np.array(frequency)
        frequency = frequency.flatten() 
        M = len(frequency)
    else:
        print("wrong data type for frequency (in getDielectricSphereFieldUnderPlaneWave)")
    
    nu = np.arange(1,N+1,1)

    [r, theta, phi] = cartToSph(sensor_location[0],sensor_location[1],sensor_location[2])

    a_n = np.ones((len(frequency), len(nu)), np.complex128)
    for c in range(0, len(nu)):
        n = nu[c]
        a_n[:,c] = (1j ** (-1*n)) * (2*n + 1) / (n * (n+1) )

    #range iterates through same numbers
    #math n not modified, index n adjusted for python zero-indexing
    aux0 = np.zeros((len(nu), 1), np.complex128); 
    aux0[0] = -1;
    aux0[1] = -3*np.cos(theta)
    for n in range(2, N):
        aux0[n] = (2*n+1)/n*np.cos(theta)*aux0[n-1] - (n+1)/n*aux0[n-2]
    
    aux1 = np.zeros((len(nu), 1), np.complex128); 
    aux1[0] = np.cos(theta);
    for n in range(2, N+1):
        aux1[n-1] = (n+1)*aux0[n-2] -n*np.cos(theta)*aux0[n-1]
    
    aux0 = np.matmul(np.ones((len(frequency),1)),np.reshape(aux0, (1, len(aux0))))
    aux1 = np.matmul(np.ones((len(frequency),1)),np.reshape(aux1, (1, len(aux1))))

    #print("aux0:\n", aux0, "\naux1:\n", aux1)
    '''
    x   = k_m*r
    print("nu:\n",nu,"\nx:\n",x)
    print(ric_bessely_derivative(nu,x,1))
    print("\nric_bessely(nu-1,x)):\n",ric_bessely(nu-1,x))
    print("\nric_bessely(0,1)):\n",ric_bessely([0],[1]))
    '''
    if r < radius:
        #not yet ready to transcribe this
        pass
    elif r > radius:
        ####calculating b_n series
        #num =  + sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselj(nu,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius)) \
        #       - sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselj_derivative(nu,k_m*radius))
        #den =  + sqrt(mu_m.*eps_s)*ones(1,N).*transpose(ric_besselj(nu,k_s*radius)).*transpose(ric_besselh_derivative(nu,2,k_m*radius))...
        #       - sqrt(mu_s.*eps_m)*ones(1,N).*transpose(ric_besselh(nu,2,k_m*radius)).*transpose(ric_besselj_derivative(nu,k_s*radius));    
        
        A =  np.matmul( np.reshape(np.sqrt(mu_s*eps_m), (M,1)), np.ones((1,N)))   
        B =  np.transpose(ric_besselj(nu,k_m*radius))
        C =  np.transpose(ric_besselj_derivative(nu,k_s*radius))

        #print(A, "\n\n", B, "\n\n", C, "\n\n")
        
        D =  np.matmul( np.reshape(np.sqrt(mu_m*eps_s), (M,1)), np.ones((1,N)))
        E =  np.transpose(ric_besselj(nu,k_s*radius))
        F =  np.transpose(ric_besselj_derivative(nu,k_m*radius))
        num =   A*B*C - D*E*F

        G = np.transpose(ric_besselh_derivative(nu, k_m*radius, 2,1))
        H = np.transpose(ric_besselh(nu,k_m*radius,2));
        den  = D*E*G - A*H*C;
        
        b_n = (num/den)*a_n
        ####calculating b_n series
        num = A*E*F - D*B*C        
        den = D*H*C - A*E*G
        c_n = (num/den)*a_n

        #print(k_m*radius)
        #print("b_n : \n", b_n, "\n\n", "c_n : \n", c_n)

        #cleaning the b_n, c_n matrices to remove inf, nan 
        # that come after values become very small
        for i in range(1, M):
            num_zeros = 0
            for j in range(0,N):
                if (abs( b_n[i,j]) < 1e-300 ):
                    num_zeros += 1
                if (num_zeros > 4):
                    b_n[i, j:] = 0
                    num_zeros = 0
                    break
        for i in range(1, M):
            num_zeros = 0
            for j in range(0,N):
                if (abs( c_n[i,j]) < 1e-300 ):
                    num_zeros += 1
                if (num_zeros > 4):
                    b_n[i, j:] = 0
                    num_zeros = 0
                    break
        
        #alpha     = (   transpose (ric_besselh_derivative (nu,2,x,2) )   +transpose(ric_besselh(nu,2,x))  )...
        #    .*transpose(  assoc_legendre(nu,1,cos(theta))   *ones(1,nFreq)  );
        x   = k_m*r
        alpha00 = np.transpose(ric_besselh_derivative(nu,x,2,2));   
        alpha01 = np.transpose(ric_besselh(nu,x,2));
        alpha10 = np.array(get_legendre(nu,1, np.cos(theta)))
        alpha11 = np.transpose( np.matmul( np.reshape(alpha10, (N,1)), np.ones((1,M)) ) )
        alpha = (alpha00 + alpha01) * alpha11

        E_r = -1j * np.cos(phi) * np.sum( (b_n * alpha) , 1)
        H_r = -1j * np.sin(phi) * np.sum( (c_n * alpha), 1) / eta_m

        alpha = np.transpose(ric_besselh_derivative(nu,x,2))* aux1
        beta = np.transpose(ric_besselh(nu,x,2)) * aux0
        summation1 = 1j*b_n*alpha - c_n*beta
        summation2 = 1j*c_n*alpha - b_n*beta
        E_theta = (np.cos(phi)/ x) * np.sum(summation1,1)
        H_theta = (np.sin(phi)/x ) * np.sum(summation2,1) / eta_m
        
        #print(ric_besselh_derivative(nu,x,2))
        #print(ric_besselj_derivative(nu,x,1))
        #print(ric_bessely_derivative(nu,x,1))
        #print("\n")
        #print(ric_bessely(nu,x,1))
        #print("theta\nalpha:\n", alpha, "\nbeta:\n", beta)

        alpha = np.transpose(ric_besselh_derivative(nu,x,2)) * aux0
        beta = np.transpose(ric_besselh(nu,x,2)) * aux1;
        summation1 = 1j*b_n*alpha - c_n*beta
        summation2 = 1j*c_n*alpha - b_n*beta
        E_phi = (np.sin(phi)/x) * np.sum(summation1,1)
        H_phi = (-1* np.cos(phi)/x) * np.sum(summation2,1) / eta_m

        #print("phi\nalpha:\n", alpha, "\nbeta:\n", beta)

        #print("E_r:\n", E_r, "\nE_theta:\n", E_theta, "\nE_phi:\n", E_phi)
        #print("H_r:\n", H_r, "\nH_theta:\n", H_theta, "\nH_phi:\n", H_phi)
    
    return [E_r, E_theta, E_phi, H_r, H_theta, H_phi]


if __name__ == '__main__':
    getDielectricSphereFieldUnderPlaneWave(0.5, DM(2.56,0.5), DM(1,0), [0,0,100], [1e6, 1e7, 1e8])
    