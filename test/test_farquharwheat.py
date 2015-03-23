# -*- coding: latin-1 -*-
"""
    test_farquhar_wheat
    ~~~~~~~~~~~~~~~~~~~
    
    Test the Farquhar-Wheat model.

    You must first install :mod:`farquharwheat` (and add it to your PYTHONPATH) 
    before running this script with the command `python`. 

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

from farquharwheat import model

def assert_close(actual, desired, tolerance=0.01):
    assert abs(actual - desired) < tolerance * (abs(desired) + 1)

def test_calculate_An():
    
    organ_width = 0.018
    organ_height = 0.6
    PAR = 24.64
    air_temperature = 15
    ambient_CO2 = 360
    humidity = 0.96
    Wind = 3.032
    organ_name = 'Lamina'
    
    actual_An, actual_Tr, _, _ = model.PhotosynthesisModel.calculate_An(organ_width, organ_height, 
        PAR, air_temperature, ambient_CO2, humidity, Wind, organ_name)
    
    desired_An = 1.22
    desired_Tr = 1.65e-6
    
    assert_close(actual_An, desired_An, tolerance=1e-3)
    assert_close(actual_Tr, desired_Tr, tolerance=1e-8)
    
    
if __name__ == '__main__':
    test_calculate_An()


