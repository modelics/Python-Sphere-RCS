import numpy as np
from DielectricMaterial import *
from src import *

class TestCase:
    '''
    '''

    def __init__(self, radius, sphere_material, background_material):
        self.radius = radius
        self.sphere_material = sphere_material
        self.background_material = background_material

class TestParameters:
    '''
    '''

    def __init__(self, sensor_location, frequency):
        self.sensor_location = sensor_location
        self.frequency = frequency

