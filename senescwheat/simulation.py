# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    senescwheat.simulation
    ~~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`senescwheat.simulation` is the front-end to run the Senesc-Wheat :mod:`model <senescwheat.model>`.
    
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


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """
    
    def __init__(self, time_step=1):
        
        #: the time step of the simulation
        self.time_step = time_step
        
        #: The inputs of Senesc-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries: 
        #:     {'roots': {(plant_index, axis_id): {roots_input_name: roots_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_id, metamer_index, organ_label, element_label): {element_input_name: element_input_value, ...}, ...}} 
        self.inputs = {}
        
        #: The outputs of Senesc-Wheat.
        #: 
        #: `outputs` is a dictionary of dictionaries: 
        #:     {'roots': {(plant_index, axis_id): {roots_output_name: roots_output_value, ...}, ...},
        #:      'elements': {(plant_index, axis_id, metamer_index, organ_label, element_label): {element_output_name: element_output_value, ...}, ...}} 
        self.outputs = {}
    
    
    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`. 
        
        :Parameters:
        
            - `inputs` (:class:`dict`) - The inputs by element.
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.  
        """
        self.inputs.clear()
        self.inputs.update(inputs)


    def run(self):
        """
        Compute Senesc-Wheat outputs from :attr:`inputs`, and update :attr:`outputs`.
        """
        self.outputs.update({inputs_type: {} for inputs_type in self.inputs.iterkeys()})
        
        all_roots_inputs = self.inputs['roots']
        all_roots_outputs = self.outputs['roots']
        for roots_inputs_id, roots_inputs_dict in all_roots_inputs.iteritems():
            mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth = model.SenescenceModel.calculate_roots_mstruct_growth(roots_inputs_dict['sucrose'], roots_inputs_dict['amino_acids'], roots_inputs_dict['mstruct'], 3600*self.time_step) #3600 is temporary, will be moved
            mstruct_death, Nstruct_death = model.SenescenceModel.calculate_roots_senescence(roots_inputs_dict['mstruct'], roots_inputs_dict['Nstruct'], 3600*self.time_step)
            delta_mstruct, delta_Nstruct, relative_delta_mstruct = model.SenescenceModel.calculate_delta_mstruct_roots(mstruct_growth, Nstruct_growth, mstruct_death, Nstruct_death, roots_inputs_dict['mstruct'])
            loss_cytokinines = model.SenescenceModel.calculate_remobilisation(roots_inputs_dict['cytokinines'], relative_delta_mstruct)
            all_roots_outputs[roots_inputs_id] = {'mstruct_C_growth': mstruct_C_growth/self.time_step, 
                                                  'Nstruct_N_growth': Nstruct_N_growth/self.time_step,
                                                  'mstruct': roots_inputs_dict['mstruct'] + delta_mstruct,
                                                  'Nstruct': roots_inputs_dict['Nstruct'] + delta_Nstruct,
                                                  'cytokinines': roots_inputs_dict['cytokinines'] - loss_cytokinines}
            
        all_elements_inputs = self.inputs['elements']
        all_elements_outputs = self.outputs['elements']
        for elements_inputs_id, elements_inputs_dict in all_elements_inputs.iteritems():
            # Senescence
            new_green_area, relative_delta_green_area, max_proteins = model.SenescenceModel.calculate_relative_delta_green_area(elements_inputs_dict['green_area'], elements_inputs_dict['proteins'] / elements_inputs_dict['mstruct'], elements_inputs_dict['max_proteins'], 3600*self.time_step)
            new_mstruct, new_Nstruct = model.SenescenceModel.calculate_delta_mstruct_shoot(relative_delta_green_area, elements_inputs_dict['mstruct'], elements_inputs_dict['Nstruct'])
            new_SLN = model.SenescenceModel.calculate_surfacic_nitrogen(elements_inputs_dict['nitrates'], elements_inputs_dict['amino_acids'], elements_inputs_dict['proteins'], elements_inputs_dict['Nstruct'], new_green_area)
            # Remobilisation
            remob_starch = model.SenescenceModel.calculate_remobilisation(elements_inputs_dict['starch'], relative_delta_green_area)
            remob_fructan = model.SenescenceModel.calculate_remobilisation(elements_inputs_dict['fructan'], relative_delta_green_area)
            remob_proteins = model.SenescenceModel.calculate_remobilisation(elements_inputs_dict['proteins'], relative_delta_green_area)
            loss_cytokinines = model.SenescenceModel.calculate_remobilisation(elements_inputs_dict['cytokinines'], relative_delta_green_area)

            # Element death and suppression from population
            min_green_area = 1E-4 # Minimal green area below which the organ is suppressed (m²)
            if new_green_area >= min_green_area:
                all_elements_outputs[elements_inputs_id] = {'green_area': new_green_area,
                                                            'mstruct': new_mstruct,
                                                            'Nstruct': new_Nstruct,
                                                            'starch': elements_inputs_dict['starch'] - remob_starch,
                                                            'sucrose': elements_inputs_dict['sucrose'] + remob_starch + remob_fructan,
                                                            'fructan': elements_inputs_dict['fructan'] - remob_fructan,
                                                            'proteins': elements_inputs_dict['proteins'] - remob_proteins,
                                                            'amino_acids': elements_inputs_dict['amino_acids'] + remob_proteins,
                                                            'cytokinines': elements_inputs_dict['cytokinines'] - loss_cytokinines}
                                                      
        
