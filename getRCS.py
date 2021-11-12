from getDielectricSphereFieldUnderPlaneWave import *
from DielectricMaterial import *
from src import *
from TestCase import *
#add any user needed imports, such as matplotlib
import matplotlib.pyplot as plt

def RCS_vs_freq(radius, ratio, background_material, sphere_material, sensor_location, save_file=None, show_plot=1):
    '''
        Calculates the RCS vs frequency for a sphere defined by 'radius' and 'sphere_material'
        located at the origin. The incident wavelengths are defined using argument 'ratio',
        where ratio is radius / wavelength. 

        Saves the plot of RCS vs freqency along with the data as a text file to 'save_file'. 

    '''
    wavelength = radius / ratio
    frequency = background_material.getPhaseVelocity(3e8 / wavelength)  / wavelength

    [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = \
        getDielectricSphereFieldUnderPlaneWave(radius, sphere_material, background_material, sensor_location, frequency)
    E = (np.stack((E_r,E_theta,E_phi), axis=0))
    mono_RCS = 4*np.pi* ( norm(sensor_location)**2 ) * np.sum( (E * np.conj(E)) , 0)

    #plotting and saving plot
    if show_plot:
        plotOneMonoRCS(radius, sphere_material, background_material, mono_RCS, frequency=frequency, savefile = save_file)
    
    #writing mono RCS data to text file, if filename is given
    if save_file:
        if not (save_file.endswith(".png")):
            save_file += ".png"
        saveMonoRCSData(save_file, mono_RCS, frequency, sphere_material, radius)

    return (frequency, mono_RCS)

def Bistatic_RCS(radius, frequency, background_material, sphere_material, distance, phi, save_filename =None, show_plot = 1):
    '''
        Calculates the bistatic RCS of a spherical object at the origin, defined by 
        sphere (DielectricMaterial) and radius, at a certain frequency.
        
        Inputs:
          radius:               radius of sphere (float, meters)
          frequency:            frequency of incident wave (float, Hz)
          background_material:  background medium (class DielectricMaterial)
          sphere_material:      sphere material (class DielectricMaterial)
          distance:             distance of observer from origin (float, meters)
          phi:                  angle of azimuth. 
          save_filename:        name of file to save the plot and data (string)

        Outputs:
          plot of bistatic RCS and the corresponding data.

        Notes:
          Assumes incident wave is travelling in the +z direction. 
          Theta = pi is the radial direction to the source point. 
          Plot saved as a .png file, data saved to .txt file.
    '''
    nang = 100
    theta = np.linspace(0,np.pi,nang)
    #distance = 2000
    #phi = 0

    bi_RCS = np.zeros((nang,), np.complex128)
    for k in range(0,nang):
        sensor_location = sphToCart(distance, theta[k], phi)
        if (type(frequency) == int or type(frequency) == float):
            frequency = np.array([frequency])

        [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = \
            getDielectricSphereFieldUnderPlaneWave(radius, sphere_material, background_material, sensor_location, frequency)
        E = (np.stack((E_r,E_theta,E_phi), axis=0))
        bi_RCS[k] = 4*np.pi* ( norm(sensor_location)**2 ) * np.sum( (E * np.conj(E)))

    #plotting and saving plot
    if save_filename:
        save_file = save_filename + ".png"
    if show_plot:
        plotBiRCS(radius, sphere_material, frequency, bi_RCS, theta, savefile = save_filename)

    #writing mono RCS data to text file, if filename is given
    if save_filename:
        data_file = save_filename + ".txt"
        saveBiRCSData(data_file, bi_RCS, theta, frequency, sphere_material)

    return (theta,bi_RCS)

def plotOneMonoRCS(radius, sphere, background, mono_RCS, *args, **kwargs):
    '''
        Plots the monostatic RCS for a single test sphere. 
        Option of saving plot to a file if specified.
        Can select different x-axis plotting presets:
            frequency (logarithmic)
            wavelength (logarigthmic)
            ratio (normal)
        
        example usage:
        sphere = DielectricMaterial(2.56,0.1, name="silicon")
        vacuum = DielectricMaterial(1,0) #background_material
        plotOneMonoRCS(radius, sphere, vacuum, my_mono_RCS_data, frequency = my_frequency_values)
        plotOneMonoRCS(radius, sphere, vacuum, my_mono_RCS_data, ratio = my_ratio, savefile = 'figure1.png')
    '''
    frequency = kwargs.get('frequency', "N/A")
    ratio = kwargs.get('ratio', "N/A")
    wavelength = kwargs.get('wavelength', "N/A")
    savefile = kwargs.get('savefile', '')
    
    if (frequency != "N/A"):
	    xseries = frequency
    elif (ratio != "N/A"):
	    xseries = ratio
    elif (wavelength != "N/A"):
	    xseries = wavelength
    else:
        print("wrong input (in plotOneMonoRCS")
        return
    
    
    plt.grid(True, which="both", ls="--")
    plt.ylabel(r'Mono-Static RCS ($m^2$)')

    if (frequency != "N/A"):
        plt.loglog(xseries, mono_RCS)
        plt.xlabel("Frequency (Hz)")
    elif (ratio != "N/A"):
        plt.semilogy(xseries, mono_RCS)
        plt.xlabel("Sphere radius in wavelengths")
    elif (wavelength != "N/A"):
        plt.loglog(xseries, mono_RCS)
        plt.xlabel("Wavelength (m)")
    

    title_str = ""

    if (sphere.sigma_e == 0):
        material = "Perfect Dielectric"
    elif (sphere.sigma_e <= 1e4):
        material = "Lossy Dielectric"
    else:
        material = "Conductor"
    
    if sphere.name:
        material = sphere.name + " Sphere"
        title_str += material
    else:
        material_params = r'$\epsilon_r$ = ' + str(round(sphere.epsilon_r,2)) + \
            r', $\mu_r$ = ' + str(round(sphere.mu_r,2)) + \
            r', $\sigma$ = ' + "{0:.2f}".format(sphere.sigma_e) + " S/m" + \
            ", radius = " + str(round(radius,2)) + " m" 
        title_str += material + " (" + material_params + ")"
    
    if (background and background.name):
        title_str += " in " + background.name
    
    plt.title(title_str)

    if (savefile):
        if not (savefile.endswith(".png")):
            savefile += ".png"
        plt.savefig(savefile, figsize=(12,6), dpi=80)  
    plt.show()

def plotBiRCS(radius, sphere, frequency, bi_RCS, theta, savefile =None):
    '''
        plots the bistatic RCS for a spherical object, defined by 
        sphere (DielectricMaterial) and radius, at a certain frequency.
        Saving the plot as an option. 
    '''
    fig, ax = plt.subplots()
    plt.semilogy(theta, bi_RCS)
    plt.grid(True, which="both", ls="--")
    
    plt.ylabel(r'Bi-Static RCS ($m^2$)')
    plt.xlabel(r'Angle $\theta$ (rad)')

    if (sphere.sigma_e == 0):
        material = "Perfect Dielectric"
    elif (sphere.sigma_e <= 1e4):
        material = "Lossy Dielectric"
    else:
        material = "Conductor"
    
    plt.title( "Bi-Static RCS for " + material + " Sphere at " + \
            "{0:.2e}".format(frequency[0]) + " Hz \n(" + \
            r' $\epsilon_r$ = ' + str(round(sphere.epsilon_r,2)) + \
            r', $\mu_r$ = ' + str(round(sphere.mu_r,2)) + \
            r', $\sigma$ = ' + '{:.2e}'.format(sphere.sigma_e) + " S/m" + \
            ", radius = " + str(round(radius,2))+ " m)"    )
    
    if (savefile):
        plt.savefig(savefile, figsize=(8,6)) 
    plt.show()

def saveMonoRCSData(savefile, mono_RCS, frequency, sphere, radius):
    '''
        Writes monostatic RCS data to a text file formatted as following.
        All data is delimited with tabs or newlines

        eps_r   ####     mu_r    ####    sigma   ####
        frequency(Hz)   RCS(m^2)
        #####   #####
        #####   #####
    '''
    if (not savefile.endswith(".txt")):
        savefile += ".txt"
    data_file = open(savefile, "w")
    header_line = "eps_r\t" + str(round(sphere.epsilon_r,2))+ \
                    "\tmu_r\t" + str(round(sphere.mu_r,2)) + \
                    "\tsigma\t" + "{0:.2e}".format(sphere.sigma_e) + \
                    "\tradius\t" + str(round(radius,2)) + "\n"
    data_file.write(header_line)
    
    column_headers = "frequency(Hz)\tRCS(m^2)\n"
    data_file.write(column_headers)

    mono_RCS = mono_RCS.flatten()
    frequency = frequency.flatten()
    n = min(len(frequency), len(mono_RCS))
    
    for i in range(0,n):
        line = "{:.9e}".format(frequency[i]) + "\t" + "{:.9e}".format(np.real(mono_RCS[i])) + "\n"
        data_file.write(line)

    data_file.close()

def saveBiRCSData(savefile, bi_RCS, theta, frequency, sphere):
    '''
        Writes bistatic RCS data to a text file formatted as following.
        All data is delimited with tabs or newlines

        freq(Hz)   ####   eps_r   ####     mu_r    ####    sigma   ####
        theta(rad)   RCS(m^2)
        #####   #####
        #####   #####
    '''
    if (not savefile.endswith(".txt")):
        savefile += ".txt"

    if (type(frequency) == list or type(frequency) == np.ndarray):
        frequency = np.array(frequency).flatten()
        frequency = float(frequency[0])


    data_file = open(savefile, "w")
    header_line =   "freq(Hz)\t" + "{0:.2e}".format(frequency) + \
                    "\teps_r\t" + str(round(sphere.epsilon_r,2))+ \
                    "\tmu_r\t" + str(round(sphere.mu_r,2)) + \
                    "\tsigma\t" + "{0:.2e}".format(sphere.sigma_e) + "\n"
    data_file.write(header_line)
    
    column_headers = "theta(rad)\tRCS(m^2)\n"
    data_file.write(column_headers)

    bi_RCS = np.reshape(bi_RCS, (bi_RCS.size,1)).flatten()
    theta = np.reshape(theta, (theta.size,1)).flatten()
    n = min(theta.size, bi_RCS.size)
    
    for i in range(0,n):
        line = "{:.9e}".format(theta[i]) + "\t" + "{:.9e}".format(np.real(bi_RCS[i])) + "\n"
        data_file.write(line)

    data_file.close()

def convertToRatio(radius, background_material, *args, **kwargs):
    '''
        given a sphere's radius and background_material 
        convert from frequency or wavelength to ratio
        return the ratio. 

        Inputs:
            radius: scalar
            background_material: of DIelectricMaterial class
            frequency: numpy array
            wavelength: numpy array

        Example: 
        input_freq = np.logspace(1,3,100)
        vacuum = DielectricMaterial(1,0)
        ratio = convertToRatio(1, vacuum, frequency = input_freq )
    '''
    c = 299792458 #speed of light
    freq = kwargs.get('frequency', 'None')
    lambd = kwargs.get('wavelength', 'None')
    rat = kwargs.get('ratio', 'None')
    if (type(rat) is not str):
        return rat
    if (type(freq) is not str):
        return freq*(radius / (background_material.getPhaseVelocity(freq)) )
    if (type(lambd) is not str):
        return radius / lambd

def Compare_RCS_vs_freq(test_cases, test_parameters, save_file = None):
    '''
        Plots several different monostatic RCS vs frequency series on the same figure.
        
        Inputs:
            test_cases: list of TestCase objects which define 
                sphere and background material, sphere radius
            test_parameters: an object of TestParameters class
                which contains information specific to the test
                such as frequency and sensor position
            save_file: filename to which to save the plot (optional)        
    '''
    fig, ax = plt.subplots()
    legend_entries = []
    for case in test_cases:
        #calcilating mono RCS as usual
        [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = \
            getDielectricSphereFieldUnderPlaneWave(case.radius, case.sphere_material, case.background_material, \
                                                    test_parameters.sensor_location, test_parameters.frequency)
        E = (np.stack((E_r,E_theta,E_phi), axis=0))
        mono_RCS = 4*np.pi* ( norm(test_parameters.sensor_location)**2 ) * np.sum( (E * np.conj(E)) , 0)
        
        #print(mono_RCS)
        #print(test_parameters.frequency)
        #plotting as usual
        plt.loglog(test_parameters.frequency, mono_RCS)
        
        series_name = "Sphere"
        if case.sphere_material.name:
            series_name = case.sphere_material.name + " " + series_name
        else:
            series_name += r'($\epsilon_r$ = ' + str(round(case.sphere_material.epsilon_r,2)) + \
             r', $\sigma$ = ' + "{0:.2f}".format(case.sphere_material.sigma_e) + " S/m)"
        if case.background_material.name:
            series_name += " in " + case.background_material.name
 
        legend_entries.append(series_name) 
    
    plt.grid(True, which="both", ls="--")
    
    plt.ylabel(r'Mono-Static RCS ($m^2$)')
    plt.xlabel("Frequency (Hz)")
    plt.legend(legend_entries, loc='best')
    plt.title("Monostatic RCS Comparison for Different Materias")

    if (save_file):
        save_filename = save_file + ".png"
        plt.savefig(save_filename, figsize=(8,6)) 
    plt.show()

def Compare_Bistatic_RCS(test_cases, test_parameters, save_file = None):
    '''
        Plots the RCS versus angle from theta = 0 to pi
        Inputs: 
          test_cases:       array of TestCase objects, length of array >=1
                            *The various spheres (with different radius, 
                            material for example) are specified within 
                            test_cases (class TestCase)
        
          test_parameters:  a single TestParameters object.
                            *The frequency, radius, and angle of azimuth phi 
                            are specified within test_parameters 
                            (class TestParameters)
          save_file:        a string denoting filename of resulting plot.
                            *Plot will be automatically saved in ".png" format
                            *If save_file not given, won't save (it is optional)
        Notes:
        * If one TestCase is given and test_parameter has only one frequency,
            the code returns the same plot as plotBiRCS()
        * If multiple TestCase objects are given within test_cases and 
          test_parameters has one frequency,
            the code returns bistatic RCS of multiple spheres at one frequency
        * If one TestCase is given and test_parameter has multiple frequencies, 
            the code performs comparison for same sphere at various frequencies 

        Example Usage:
            vacuum = DielectricMaterial(1,0)
            case1 = TestCase(0.5, DielectricMaterial(2.56,0), vacuum)
            case2 = TestCase(0.5, DielectricMaterial(2.56,0.1, name = "Silicon"), vacuum)

            distance = 2000
            phi = np.pi/2
            sensor_location = cartToSph(distance, theta=0, phi)
            frequency = 1e9
            param1 = TestParameters(sensor_location, frequency)

            Compare_RCS_vs_freq([case1,case2], param1, "two_material_comparison_example")
    '''
    #checking the frequency input
    frequency = test_parameters.frequency
    if (type(frequency) == int or type(frequency) == float):
        frequency = np.array([frequency])
    if (type(frequency) == list or type(frequency) == np.ndarray):
        frequency = np.array(frequency).flatten()
    
    #just a single bistatic RCS plot
    if (len(test_cases) == 1 and frequency.size == 1 ):
        print("Compare_Bistatic_RCS cases 1, frequency 1")
        case = test_cases[0]
        param = test_parameters[0]
        radius = case.radius
        background = case.background_material
        sphere = case.sphere_material
        sphere_location = param.sensor_location
        frequency = param.frequency

        [distance, theta, phi] = sphToCart(sensor_location)
        #generate RCS vs theta data
        (theta, bi_RCS) = Bistatic_RCS(radius, frequency, background, sphere, distance, phi, save_filename =None, show_plot=0)
        # plot the data
        plotBiRCS(radius, sphere, frequency, bi_RCS, theta, savefile =save_file)
    
    #multiple spheres, one frequency
    if (len(test_cases) >= 1 and frequency.size == 1):
        print("Compare_Bistatic_RCS cases many, frequency 1")
        fig, ax = plt.subplots()
        legend_entries = []

        for case in test_cases:
            radius = case.radius
            background = case.background_material
            sphere = case.sphere_material

            x = test_parameters.sensor_location[0]
            y = test_parameters.sensor_location[1]
            z = test_parameters.sensor_location[2]
            [distance, phi, theta] = cartToSph(x,y,z)

            (theta, bi_RCS) = Bistatic_RCS(radius, frequency, background, sphere, distance, phi, show_plot=0)               
            plt.semilogy(theta, bi_RCS)
            
            series_name = "Sphere"
            if case.sphere_material.name:
                series_name = case.sphere_material.name + " " + series_name
            else:
                series_name += r'($\epsilon_r$ = ' + str(round(case.sphere_material.epsilon_r,2)) + \
                r', $\sigma$ = ' + "{0:.2f}".format(case.sphere_material.sigma_e) + " S/m)"
            if case.background_material.name:
                series_name += " in " + case.background_material.name
    
            legend_entries.append(series_name) 
        
        plt.ylabel(r'Bi-Static RCS ($m^2$)')
        plt.xlabel(r'Angle $\theta$ (rad)')
        plt.grid(True, which="both", ls="--")
        
        plt.legend(legend_entries, loc='best')
        plt.title("Bistatic RCS Comparison for Different Materias at " + "{0:.2e}".format(frequency) + " Hz")

        if (save_file):
            save_filename = save_file + ".png"
            plt.savefig(save_filename, figsize=(8,6)) 
        plt.show()
    
    #one sphere, multiple frequency
    elif(len(test_cases) == 1 and frequency.size > 1):
        print("Compare_Bistatic_RCS cases 1, frequency many")
        case = test_cases[0]
        radius = case.radius
        background = case.background_material
        sphere = case.sphere_material

        x = test_parameters.sensor_location[0]
        y = test_parameters.sensor_location[1]
        z = test_parameters.sensor_location[2]
        [distance, phi, theta] = cartToSph(x,y,z)

        fig, ax = plt.subplots()
        legend_entries = []

        for k in range(frequency.size):
            #print("frequency: " , frequency)
            #print("frequency[k]: ", frequency[k], type(frequency[k]))
            (theta, bi_RCS) = Bistatic_RCS(radius, frequency[k:k+1], background, sphere, distance, phi, show_plot=0)               
            plt.semilogy(theta, bi_RCS)
            
            series_name = "Sphere"
            if case.sphere_material.name:
                series_name = case.sphere_material.name + " " + series_name
            else:
                series_name += r'($\epsilon_r$ = ' + str(round(case.sphere_material.epsilon_r,2)) + \
                r', $\sigma$ = ' + "{0:.2f}".format(case.sphere_material.sigma_e) + " S/m)"
            if case.background_material.name:
                series_name += " in " + case.background_material.name
            series_name += " at "+ "{0:.2e}".format(frequency[k]) + " Hz"
    
            legend_entries.append(series_name) 
        
        plt.ylabel(r'Bi-Static RCS ($m^2$)')
        plt.xlabel(r'Angle $\theta$ (rad)')
        plt.grid(True, which="both", ls="--")
        
        plt.legend(legend_entries, loc='best')
        plt.title("Bistatic RCS Comparison for Different Frequencies")

        if (save_file):
            save_filename = save_file + ".png"
            plt.savefig(save_filename, figsize=(8,6)) 
        plt.show() 
        
    #extra: catching errors
    elif (len(test_cases) > 1 and frequency.size > 1):
        print("inputs not supported")
    else:
        print("error in compare_bistatic_rcs")

def plotFromFile(filenames, plt_type, save_file=''):
    '''
        Reads RCS data from filenames and displays on the same plot.

        Inputs:
            filenames:  a list of strings. strings are the file name  
                        of the RCS data with  file extension ".txt"
            plt_type:   string. if type = 'bi' or 'bistatic', returns 
                        bistatic RCS plot. otherwise, if 'mono', or
                        'monostatic', returns monostatic plot.
            save_file:  If given, saves plot to this filename as png 
        
        Note: 
            * Data files are assumed to be in the same format as this 
              code exports. See saveMonoRCSData() and saveBiRCSData()
              for details.
    '''
    plt_type = plt_type.replace(" ", "").lower()

    if (plt_type == 'mono' or plt_type == 'monostatic'):
        fig, ax = plt.subplots()
        legend_entries = []
        
        for filename in filenames:
            data_file = open(filename, "r")
            line_1 = data_file.readline()
            # line 1:
            #eps_r	2.56	mu_r	1	sigma	3.00e-02
            [eps_r, mu_r, sigma, radius] = line_1.split()[1::2]
            sphere = DielectricMaterial(float(eps_r), float(sigma), float(mu_r),0)
            radius = float(radius)

            #skip next line which is column header
            #frequency(Hz)	RCS(m^2)
            data_file.readline()

            frequency, monoRCS = [],[]
            lines = data_file.readlines()
            for line in lines:
                values = line.split()
                frequency.append(float(values[0]))
                monoRCS.append(float(values[1]))
            
            frequency = np.array(frequency).flatten()
            monoRCS = np.array(monoRCS).flatten()
            data_file.close()

            #print("frequency: \n", frequency)
            #print("\nmonoRCS: \n", monoRCS)
            #print("\n just before plot")
            
            plt.loglog(frequency, monoRCS)

            series_name = ""   
            if sphere.name:
                material = sphere.name
                series_name += material + "radius = " + str(round(radius,2)) + " m)"         
            else:
                descriptor = r'($\epsilon_r$ = ' + str(round(sphere.epsilon_r,2)) + \
                        r', $\sigma$ = ' + "{0:.2f}".format(sphere.sigma_e) + " S/m" + \
                        ", radius = " + str(round(radius,2)) + " m)"
                series_name +=  descriptor
            
            legend_entries.append(series_name) 
        
        plt.grid(True, which="both", ls="--")
        plt.ylabel(r'Mono-Static RCS ($m^2$)')
        plt.xlabel("Frequency (Hz)")
        plt.legend(legend_entries, loc='best')
        plt.title("Monostatic RCS Comparison for Different Materias")

        if save_file:
            save_file += ".png"
            plt.savefig(save_file, figsize=(8,6))
        plt.show()

    elif (plt_type == 'bi' or plt_type == 'bistatic'):
        fig, ax = plt.subplots()
        legend_entries = []
        
        for filename in filenames:
            data_file = open(filename, "r")
            line_1 = data_file.readline()
            # line 1:
            #eps_r	2.56	mu_r	1	sigma	3.00e-02
            [frequency, eps_r, mu_r, sigma] = line_1.split()[1::2]
            sphere = DielectricMaterial(float(eps_r), float(sigma), float(mu_r),0)
            frequency = float(frequency)

            #skip next line which is column header
            #frequency(Hz)	RCS(m^2)
            data_file.readline()

            theta, biRCS = [],[]
            lines = data_file.readlines()
            for line in lines:
                values = line.split()
                theta.append(float(values[0]))
                biRCS.append(float(values[1]))
            
            theta = np.array(theta).flatten()
            biRCS = np.array(biRCS).flatten()
            data_file.close()

            #print("theta: \n", theta)
            #print("\nbiRCS: \n", biRCS)
            #print("\n just before plot")
            
            plt.semilogy(theta, biRCS)

            series_name = ""   
            if sphere.name:
                material = sphere.name
                series_name += material #+ "radius = " + str(round(radius,2)) + " m)"         
            else:
                descriptor = r'Sphere ($\epsilon_r$ = ' + str(round(sphere.epsilon_r,2)) + \
                        r', $\sigma$ = ' + "{0:.2f}".format(sphere.sigma_e) + " S/m)" +\
                        " at " + "{0:.2e}".format(frequency) + "Hz"
                        #", radius = " + str(round(radius,2)) + " m)"
                series_name +=  descriptor
            
            legend_entries.append(series_name) 
        
        plt.grid(True, which="both", ls="--")
        plt.ylabel(r'Bi-Static RCS ($m^2$)')
        plt.xlabel(r'Angle $\theta$ (rad)')
        plt.legend(legend_entries, loc='best')
        plt.title("Bistatic RCS Comparison from Data Files")

        if save_file:
            save_file += ".png"
            plt.savefig(save_file, figsize=(8,6))
        plt.show()
    else:
        print("incorrect plt_type input to plotFromFile()")


if __name__ == '__main__':
    '''
    radius = 0.5 #meters
    ratio = np.arange(0.01,1.61,0.01)
    wavelength = radius / ratio

    background = DielectricMaterial(1,0)
    frequency = background.getPhaseVelocity(3e8 / wavelength)  / wavelength
    
    sensor_location = [0,0,-1000]
    sphere = DielectricMaterial(2.56, 0.0)
    sphere = DielectricMaterial(1e8,0,1e-8,0)
   
    [E_r, E_theta, E_phi, H_r, H_theta, H_phi] = \
        getDielectricSphereFieldUnderPlaneWave(radius, sphere, background, sensor_location, frequency)
    E = (np.stack((E_r,E_theta,E_phi), axis=0))
    mono_RCS = 4*np.pi* ( norm(sensor_location)**2 ) * np.sum( (E * np.conj(E)) , 0)
    '''


    #print(mono_RCS)

    #plotOneMonoRCS(radius, sphere, background, mono_RCS, ratio = ratio, savefile="PEC_ratio")

    #(theta, bi_RCS) = Bistatic_RCS(radius, 1e9, background, sphere, 2000, 0, show_plot=0)
    #saveBiRCSData("bistatic_perfect_dielectric_example", bi_RCS, theta, 1e9, sphere)

    #print(convertToRatio(radius, background, wavelength = wavelength))

    #print(background.name)

    #testing sigma sweep, or material comparison
    '''
        case1 = TestCase(0.5, DielectricMaterial(2.56,0), DielectricMaterial(1,0))
        case2 = TestCase(0.5, DielectricMaterial(2.56,0.1, name = "Silicon"), DielectricMaterial(1,0))
        param1 = TestParameters([0,0,-2000], np.logspace(7,9,100))
        Compare_RCS_vs_freq([case1,case2], param1, "example_2_material_comparison")
    '''
    
    #testing PEC sphere
    #(freq, mono_RCS) = RCS_vs_freq(radius = 0.5, ratio = np.arange(0.01,1.61,0.01), background_material = DielectricMaterial(1,0), sphere_material = DielectricMaterial(1e8,0,1e-8,0), sensor_location = [0,0,-2000], save_file = 'PEC', show_plot = 1)

    #testing plot with ratio on x axis
    '''
        radius = 0.5 #meters
        ratio = np.arange(0.01,1.61,0.01)
        background = DielectricMaterial(1,0)
        sphere_material = DielectricMaterial(2.56,0.003)
        sensor_location = [0,0,-2000]

        (freq, mono_RCS) = RCS_vs_freq(radius, ratio, background, sphere_material, \
                                sensor_location , save_file = 'PEC', show_plot = 0)
        plotOneMonoRCS(radius, sphere_material, background, mono_RCS, ratio = ratio, \
                        savefile = "lossy_dielectric_mono_rcs")
    '''

    #testing perfect dieelctric sphere
    #RCS_vs_freq(radius = 0.5, ratio = np.arange(0.01,1.61,0.01), background_material = DielectricMaterial(1,0), sphere_material = DielectricMaterial(2.56,0), sensor_location = [0,0,-2000], save_file = 'perfect_dielectric')

    #testing lossy dielectric sphere
    #RCS_vs_freq(radius = 0.5, ratio = np.arange(0.01,1.61,0.01), background_material = DielectricMaterial(1,0), sphere_material = DielectricMaterial(2.56,0.03), sensor_location = [0,0,-2000], save_file = 'lossy_dielectric')
    
    #stress-testing getNMax function: full sigma sweep
    '''
        #auto-generating test cases from list of conductivities
        vacuum = DielectricMaterial(1,0)
        radius = 0.5
        conductivities = [0, 1e-5, 1e-3, 1e-2, 1e0, 1e3, 1e6]
        names  = ["0", "1e-5", "1e-3", "1e-2", "1e0", "1e3", "1e6"]
        #conductivities = [0,1e0,1e6]
        #names = ["0", "1e0", "1e6"]
        eps_r = 2.56
        mu_r = 1
        cases = []
        for i in range(0, len(conductivities)):
            sphere = DielectricMaterial(eps_r, conductivities[i], mu_r, name = names[i])
            test_case = TestCase(radius, sphere, vacuum)
            cases.append(test_case)
        
        #formalizing test parameters
        sensor_location = [0,0,-2000]
        frequency = np.logspace(7,9,1000)
        param1 = TestParameters(sensor_location, frequency)

        Compare_RCS_vs_freq(cases, param1, save_file = "sigma_sweep")
    '''

    #testing Compare_Bistatic_RCS
    '''
        vacuum = DielectricMaterial(1,0)
        radius = 0.5
        sphere1 = DielectricMaterial(2.56,0, 1,0, name = "Loss-less")
        sphere2 = DielectricMaterial(2.56,1, 1,0, name = "Lossy")
        test_cases = [TestCase(radius, sphere1, vacuum), TestCase(radius, sphere2, vacuum)]

        sensor_location = [0,0,-2000]
        frequency = 1e9
        test_parameters = TestParameters(sensor_location, frequency)
        
        #bistatic RCS for two spheres, one frequency
        #   Compare_Bistatic_RCS(test_cases, test_parameters, save_file = "compare_bistatic_materials")

        #bistatic RCS for one sphere, two frequencies
        vacuum = DielectricMaterial(1,0)
        radius = 0.5
        sphere3 = DielectricMaterial(2.56,3.3, 1,0, name = "Lossy")
        test_cases = [TestCase(radius, sphere3, vacuum)]

        sensor_location = [0,0,-2000]
        frequency = [1e9, 3e9, 5e9]
        test_parameters = TestParameters(sensor_location, frequency)

        Compare_Bistatic_RCS(test_cases, test_parameters, save_file = "compare_bistatic_frequencies")
    '''

        
    #DEFAULT PLOTTING CODE
    '''
        plt.loglog(frequency,mono_RCS)
        plt.grid(True, which="both", ls="--")
        
        plt.ylabel(r'Mono-Static RCS ($m^2$)')
        plt.xlabel("Frequency (Hz)")

        if (sphere.sigma_e == 0):
            material = "Perfect Dielectric"
        elif (sphere.sigma_e <= 1e4):
            material = "Lossy Dielectric"
        else:
            material = "Conductor"
        
        plt.title(material + r' $\epsilon_r$ = ' + str(round(sphere.epsilon_r,2)) + \
                r', $\mu_r$ = ' + str(round(sphere.mu_r,2)) + \
                r', $\sigma$ = ' + "{0:.2f}".format(sphere.sigma_e) + " S/m" + \
                ", radius = " + str(round(radius,2)) + " m" )
        plt.show()
    '''

    #testing plot from file
    #plotFromFile(['lossy_dielectric.txt', 'perfect_dielectric.txt'], 'monostatic', save_file='plot_from_file_example')
    #plotFromFile(['bistatic_lossy_dielectric_example.txt', "bistatic_perfect_dielectric_example.txt"], 'bi', save_file='plot_from_file_example_2')

    #validating with MoM Solver Results
    radius = 0.5
    frequency = 1000000000
    background = DielectricMaterial(1,0)
    sphere = DielectricMaterial(2,1000)
    distance = 2000
    phi = 0

    Bistatic_RCS(radius, frequency, background, sphere, distance, phi, save_filename ="calculations", show_plot = 0)
    plotFromFile(['sigma1e3_eaefie_bRCS_f_1000000000_Hz.txt', 'calculations.txt'], 'bi')
    #old debugging inputs
    '''
    radius = 0.5;
    sphere = DielectricMaterial(2.56,0.5)
    background = DielectricMaterial(1,0)
    sensor_location = [0,0,100];
    frequency = [1e6, 1e7, 1e8]
    '''