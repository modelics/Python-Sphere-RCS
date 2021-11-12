import numpy as np
import math
from DielectricMaterial import DielectricMaterial as DM
from src import *
from bessel import *

def getNMax(radius, sphere, background, frequency):
    '''
        determines the appropriate number of mie terms to evaluate
        Based on the wiscombe 1980 recommendation (which was deternined through 
        convergence behavior bessel functions at high orders that were
        calculated recursivelly).

        Designed to work for single-layered (monolithic) sphere. 
    '''
    #check that frequency input is correct  
    if (type(frequency) == int or type(frequency) == float):
        frequency = np.array([frequency])
    if (type(frequency) == list or type(frequency) == np.ndarray):
        frequency = np.array(frequency).flatten()
        M = len(frequency)
    else:
        print("wrong data type for frequency (in getNMax)")
    
    
    k_m = DM.getWaveNumber(background, frequency)
    x = abs(k_m * radius)
    #print(x)

    N_m = DM.getComplexRefractiveIndex(background, frequency)
    m = DM.getComplexRefractiveIndex(sphere, frequency) / N_m #relative refractive index

    N_max = np.ones((M,))
    for k in range(0,M):
        if (0.02 <= x[k] and x[k] <= 8):
            N_stop = x[k] + 4.*x[k]**(1/3) + 1
        elif (8 < x[k] and x[k] < 4200):
            N_stop = x[k] + 4.05*x[k]**(1/3) + 2
        elif (4200 <= x[k] and x[k] <= 20000):
            N_stop = x[k] + 4.*x[k]**(1/3) + 2
        else:
            print("WaRNING: size parameter over 20,000")
            N_stop = 20000 + 4.*20000**(1/3) + 2
        
        #this is the KZHU original nmax formula (adapted for single sphere)
        #it recommends 100's of terms for real metals
        #N_max[k] = max(N_stop, abs(m[k] * x[k]) )+15

        #this is the Wiscombe-only implementation, seems to be accurate enough
        N_max[k] = N_stop
        
    return math.ceil(max(N_max))


if __name__ == "__main__":
    radius = 0.5
    sphere = DM(5,100)
    background = DM(1,0)
    frequency = np.logspace(7,9,5)
    print(getNMax(radius, sphere, background, frequency))

    print(sphere.getComplexRefractiveIndex(1e9))