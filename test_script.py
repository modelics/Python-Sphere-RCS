'''
    Purpose of document: run this in the terminal to confirm 
    that all functions in this library work as expected.
'''

from getRCS import *
import warnings

if __name__ == "__main__":
    print("------------------\nPYTHON SPHERE RCS\n------------------\n")
    print("this script will confirm that all functions work as expected.")
    print("It will display plots and save them to the \"output\" folder \n")
    confirm = input("As you see plots appearing, close the Matplotlib window to proceed. Hit Enter to Continue:")
    print("------\nBEGIN\n------\n")

    #PART 1: PLOT ONE MONOSTATIC RCS FOR LOSSY DIELECTRIC
    radius = 0.5 #meters
    ratio = np.arange(0.01,1.61,0.01)
    background = DielectricMaterial(1,0)
    sphere_material = DielectricMaterial(2.56,0.003)
    sensor_location = [0,0,-2000]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        (freq, mono_RCS) = RCS_vs_freq(radius, ratio, background, sphere_material, \
                            sensor_location , save_file = 'PEC', show_plot = 0)
        plotOneMonoRCS(radius, sphere_material, background, mono_RCS, ratio = ratio, \
                    savefile = "output/lossy_dielectric_mono_rcs")
    
    print("Finished Part 1: Plot Monostatic RCS for Losy Dielectric Sphere as lossy_dielectric_mono_rcs.png\n")

    #PART 2: PLOT ONE MONOSTATIC RCS FOR PERFECT DIELECTRIC
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        RCS_vs_freq(radius = 0.5, ratio = np.arange(0.01,1.61,0.01), \
            background_material = DielectricMaterial(1,0), sphere_material = DielectricMaterial(2.56,0), \
            sensor_location = [0,0,-2000], save_file = 'output/perfect_dielectric')
    
    print("Finished Part 2: Plot Monostatic RCS for Perfect Dielectric Sphere as perfect_dielectric.png\n")

    #PART 3A: PLOT BISTATIC RCS COMPARISON 
    #bistatic RCS for two spheres, one frequency
    vacuum = DielectricMaterial(1,0)
    radius = 0.5
    sphere1 = DielectricMaterial(2.56,0, 1,0, name = "Loss-less")
    sphere2 = DielectricMaterial(2.56,1, 1,0, name = "Lossy")
    test_cases = [TestCase(radius, sphere1, vacuum), TestCase(radius, sphere2, vacuum)]

    sensor_location = [0,0,-2000]
    frequency = 1e9
    test_parameters = TestParameters(sensor_location, frequency)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Compare_Bistatic_RCS(test_cases, test_parameters, save_file = "output/compare_bistatic_materials")
    
    #PART 3B: PLOT BISTATIC RCS COMPARISON 
    #bistatic RCS for one sphere, two frequencies
    vacuum = DielectricMaterial(1,0)
    radius = 0.5
    sphere3 = DielectricMaterial(2.56,3.3, 1,0, name = "Lossy")
    test_cases = [TestCase(radius, sphere3, vacuum)]

    sensor_location = [0,0,-2000]
    frequency = [1e9, 3e9, 5e9]
    test_parameters = TestParameters(sensor_location, frequency)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        Compare_Bistatic_RCS(test_cases, test_parameters, save_file = "output/compare_bistatic_frequencies")

    print("Finished Part 3: Compare Bistatic RCS for Perfect Dielectric Spheres, parts A and B\n")

    print("------\nFINISHED\n------\n")