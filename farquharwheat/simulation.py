# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.simulation
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.simulation` is the front-end to run the Farquhar-Wheat :mod:`model <farquharwheat.model>`.
    
    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
"""

"""
    Information about this versioned file:
        $LastChangedBy$
        $LastChangedDate$
        $LastChangedRevision$
        $URL$
        $Id$
"""

import pandas as pd

import model

class Simulation(object):
    """
    The Simulation class permits to initialize and run the model, and format the 
    outputs of the model.

    Use:
    
        * :meth:`initialize` to initialize the model,
        * :meth:`run` to run the model,
        * :meth:`format_outputs` to format the outputs of the model.

    """
    
    ELEMENTS_KEYS_NAMES = ['plant', 'axis', 'metamer', 'organ', 'element'] #: the keys which define the topology of an element.
    
    def __init__(self):
        #: The inputs by element.
        #:
        #: `inputs` is a dictionary of dictionaries: {element1_id: element1_inputs, element2_id: element2_inputs, ..., elementN_id: elementN_inputs}, where: 
        #:     * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_type, element_type),
        #:     * and elementi_inputs is a dictionary: {'elementi_input1_name': elementi_input1_value, 'elementi_input2_name': elementi_input2_value, ..., 'elementi_inputN_name': elementi_inputN_value}.
        #: 
        #: See :meth:`PhotosynthesisModel.calculate_An <farquharwheat.model.PhotosynthesisModel.calculate_An>` 
        #: for more information about the inputs.  
        self.inputs = {}
        #: The outputs by element. 
        #: 
        #: `outputs` is a dictionary of dictionaries: {element1_id: element1_outputs, element2_id: element2_outputs, ..., elementN_id: elementN_outputs}, where: 
        #:     * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_type, element_type),
        #:     * and elementi_outputs is a dictionary: {'elementi_output1_name': elementi_output1_value, 'elementi_output2_name': elementi_output2_value, ..., 'elementi_outputN_name': elementi_outputN_value}.
        #: 
        #: See :meth:`PhotosynthesisModel.calculate_An <farquharwheat.model.PhotosynthesisModel.calculate_An>` 
        #: for more information about the outputs.
        self.outputs = {}
    
    
    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`. 
        
        :Parameters:
        
            - `inputs` (:class:`dict` of :class:`dict`) - The inputs by element.
                `inputs` is a dictionary of dictionaries: {element1_id: element1_inputs, element2_id: element2_inputs, ..., elementN_id: elementN_inputs}, where:
                 
                    * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_type, element_type),
                    * and elementi_inputs is a dictionary: {'elementi_input1_name': elementi_input1_value, 'elementi_input2_name': elementi_input2_value, ..., 'elementi_inputN_name': elementi_inputN_value}.
                    
                See :meth:`PhotosynthesisModel.calculate_An <farquharwheat.model.PhotosynthesisModel.calculate_An>` 
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
            organ_type = element_inputs['organ_type']
            Na = element_inputs['Na']
            organ_width = element_inputs['organ_width']
            organ_height = element_inputs['organ_height']
            STAR = element_inputs['STAR']
            PAR = STAR * PARi
            Ag, An, Rd, Tr, Torg, gs = model.PhotosynthesisModel.calculate_An(Na, organ_width, organ_height, PAR, Ta, ambient_CO2, RH, Ur, organ_type)
            self.outputs[element_id] = {'Ag': Ag, 'An': An, 'Rd': Rd, 'Tr': Tr, 'Torg': Torg, 'gs': gs}
    
            
    def format_outputs(self): 
        """
        Format :attr:`outputs` to Pandas dataframe.
        
        :Returns:
            The outputs in a dataframe, with one row by element.
            See :meth:`PhotosynthesisModel.calculate_An <farquharwheat.model.PhotosynthesisModel.calculate_An>` 
            for more information about the outputs.  
        
        :Returns Type:
            :class:`pandas.DataFrame`
        
        """
        elements_ids_df = pd.DataFrame(self.outputs.keys(), columns=Simulation.ELEMENTS_KEYS_NAMES)
        elements_outputs_df = pd.DataFrame(self.outputs.values())
        outputs_df = pd.concat([elements_ids_df, elements_outputs_df], axis=1)
        outputs_df.sort_index(by=Simulation.ELEMENTS_KEYS_NAMES, inplace=True)
        return outputs_df
            
            
