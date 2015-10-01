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
        #: The inputs by element.
        #:
        #: `inputs` is a dictionary of dictionaries: {element1_id: element1_inputs, element2_id: element2_inputs, ..., elementN_id: elementN_inputs}, where: 
        #:     * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_label, element_type),
        #:     * and elementi_inputs is a dictionary: {'elementi_input1_name': elementi_input1_value, 'elementi_input2_name': elementi_input2_value, ..., 'elementi_inputN_name': elementi_inputN_value}.
        #: 
        #: See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
        #: for more information about the inputs.  
        self.inputs = {}
        #: The outputs by element. 
        #: 
        #: `outputs` is a dictionary of dictionaries: {element1_id: element1_outputs, element2_id: element2_outputs, ..., elementN_id: elementN_outputs}, where: 
        #:     * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_label, element_type),
        #:     * and elementi_outputs is a dictionary: {'elementi_output1_name': elementi_output1_value, 'elementi_output2_name': elementi_output2_value, ..., 'elementi_outputN_name': elementi_outputN_value}.
        #: 
        #: See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
        #: for more information about the outputs.
        self.outputs = {}
    
    
    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`. 
        
        :Parameters:
        
            - `inputs` (:class:`dict`) - The inputs by element.
               Inputs must be a dictionary of dictionaries: {element1_id: element1_inputs, element2_id: element2_inputs, ..., elementN_id: elementN_inputs}, where:
                 
                   * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_label, element_type),
                   * and elementi_inputs is a dictionary: {'elementi_input1_name': elementi_input1_value, 'elementi_input2_name': elementi_input2_value, ..., 'elementi_inputN_name': elementi_inputN_value}.
                    
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
            
            - `PARi` (:class:`float`) - the incident PAR (µmol m-2 s-1)
            
        """
        self.outputs.clear()
        for (element_id, element_inputs) in self.inputs.iteritems():
            organ_label = element_inputs['label']
            surfacic_nitrogen = element_inputs['surfacic_nitrogen']
            width = element_inputs['width']
            height = element_inputs['height']
            STAR = element_inputs['STAR']
            PAR = STAR * PARi
            Ag, An, Rd, Tr, Ts, gs = model.Model.calculate_An(surfacic_nitrogen, width, height, PAR, Ta, ambient_CO2, RH, Ur, organ_label)
            self.outputs[element_id] = {'Ag': Ag, 'An': An, 'Rd': Rd, 'Tr': Tr, 'Ts': Ts, 'gs': gs}
    
