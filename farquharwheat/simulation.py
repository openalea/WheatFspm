# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.simulation
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.simulation` is the front-end to run the Farquhar-Wheat :mod:`model <farquharwheat.model>`.
    
    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import model

class SimulationError(Exception): pass
class SimulationInputsError(SimulationError): pass

class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """
    
    def __init__(self):
        
        #: The inputs of Farquhar-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries: 
        #:     {'organs': {(plant_index, axis_label): {organs_input_name: organs_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_input_name: element_input_value, ...}, ...}}
        #: See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
        #: for more information about the inputs.  
        self.inputs = {}
        
        #: The outputs of Farquhar-Wheat.
        #: 
        #: `outputs` is a dictionary of dictionaries: 
        #:     {'organs': {(plant_index, axis_label): {organs_output_name: organs_output_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_output_name: element_output_value, ...}, ...}}
        #: See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
        #: for more information about the outputs. 
        self.outputs = {}
        
    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`. 
        
        :Parameters:
        
            - `inputs` (:class:`dict`) - The inputs by organ and element.
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.
              
            See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
               for more information about the inputs.
        """
        self.inputs.clear()
        self.inputs.update(inputs)


    def run(self, Ta, ambient_CO2, RH, Ur, PARi):
        """
        Compute Farquhar variables for each element in :attr:`inputs` and put 
        the results in :attr:`outputs`.
        
        :Parameters:
        
            - `Ta` (:class:`float`) - air temperature at t (degree Celsius)

            - `ambient_CO2` (:class:`float`) - air CO2 at t (umol mol-1)

            - `RH` (:class:`float`) - relative humidity at t (decimal fraction)

            - `Ur` (:class:`float`) - wind speed at the top of the canopy at t (m s-1)
            
            - `PARi` (:class:`float`) - the incident PAR above the canopy (µmol m-2 s-1)
            
        """
        self.outputs.update({inputs_type: {} for inputs_type in self.inputs.iterkeys()})
        organs_inputs = self.inputs['organs']
        organs_outputs = self.outputs['organs']
        elements_inputs = self.inputs['elements']
        elements_outputs = self.outputs['elements']
        for (element_id, element_inputs) in elements_inputs.iteritems():
            organ_id = tuple(element_id[:-1])
            if element_inputs['green_area'] == 0.0:
                element_outputs = dict.fromkeys(['Ag', 'An', 'Rd', 'Tr', 'Ts', 'gs'], 0.0)
            else:
                organ_label = organs_inputs[organ_id]['label']
                STAR = element_inputs['STAR'] # TODO: check whether absorbed STAR or not.
                PARa = STAR * PARi
                surfacic_nitrogen = model.Model.calculate_surfacic_nitrogen(element_inputs['nitrates'], 
                                                                            element_inputs['amino_acids'], 
                                                                            element_inputs['proteins'], 
                                                                            element_inputs['Nstruct'], 
                                                                            element_inputs['green_area'])
                Ag, An, Rd, Tr, Ts, gs = model.Model.calculate_An(surfacic_nitrogen, 
                                                                  element_inputs['width'], 
                                                                  element_inputs['height'], 
                                                                  PARa, Ta, ambient_CO2, RH, Ur, organ_label)
                element_outputs = {'Ag': Ag, 'An': An, 'Rd': Rd, 'Tr': Tr, 'Ts': Ts, 'gs': gs}
                
            organs_outputs[organ_id] = {}
            elements_outputs[element_id] = element_outputs
