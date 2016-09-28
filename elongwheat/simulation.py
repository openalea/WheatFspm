# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    elongwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`elongwheat.simulation` is the front-end to run the CN-Wheat model.

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

import numpy as np
import pandas as pd

import model
import logging
import warnings
import copy

#: the inputs needed by ElongWheat
HIDDENZONE_INPUTS = ['leaf_is_growing', 'hiddenzone_L', 'leaf_L', 'leaf_Lmax', 'leaf_Lem_prev', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'mstruct']
ORGAN_INPUTS = ['visible_length', 'is_growing', 'final_hidden_length', 'length']

#: the outputs computed by ElongWheat
HIDDENZONE_OUTPUTS = ['leaf_is_growing', 'hiddenzone_L', 'leaf_L', 'delta_leaf_L', 'leaf_Lmax', 'leaf_Lem_prev', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'mstruct']
ORGAN_OUTPUTS = ['visible_length', 'is_growing', 'final_hidden_length', 'length']

#: the inputs and outputs of ElongWheat.
HIDDENZONE_INPUTS_OUTPUTS = sorted(set(HIDDENZONE_INPUTS + HIDDENZONE_OUTPUTS))
ORGAN_INPUTS_OUTPUTS = sorted(set(ORGAN_INPUTS + ORGAN_OUTPUTS))


class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, delta_t=1):

        #: The inputs of elong-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'hiddenzone_L_calculation': {(plant_index, axis_label, metamer_index): {'previous_hiddenzone_length': previous_hiddenzone_length,
        #:                                                                       'previous_sheath_visible_length': previous_sheath_visible_length,
        #:                                                                       'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of elong-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'organs': {(plant_index, axis_label, metamer_index, organ_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'hiddenzone_L_calculation': {(plant_index, axis_label, metamer_index): {'previous_hiddenzone_length': previous_hiddenzone_length,
        #:                                                                       'previous_sheath_visible_length': previous_sheath_visible_length,
        #:                                                                       'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.outputs = {}

        #: the delta t of the simulation (in seconds)
        self.delta_t = delta_t


    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :Parameters:

            - `inputs` (:class:`dict`)
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self):
        # Copy the inputs into the output dict
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hiddenzone', 'organs', 'hiddenzone_L_calculation'])})

        # Hidden zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # organs
        all_organs_inputs = self.inputs['organs']
        all_organs_outputs = self.outputs['organs']

        # Previous sheaths
        all_hiddenzone_L_calculation_inputs = self.inputs['hiddenzone_L_calculation']

        for hiddenzone_id, hiddenzone_inputs in all_hiddenzone_inputs.iteritems():
            curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]
            if hiddenzone_id==(1,'MS',4):
                pass
            # Found previous hidden zone TODO: a améliorer
            prev_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] - 1])
            if prev_hiddenzone_id in all_hiddenzone_inputs:
                prev_leaf_emerged = all_hiddenzone_inputs[prev_hiddenzone_id]['leaf_is_emerged']
            else:
                prev_leaf_emerged = True

            # Update hidden zone length
            previous_hiddenzone_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['previous_hiddenzone_length']
            previous_sheath_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['previous_sheath_visible_length']
            previous_sheath_final_hidden_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['previous_sheath_final_hidden_length']
            curr_hiddenzone_outputs['hiddenzone_L'] = model.calculate_hiddenzone_length(previous_hiddenzone_L, previous_sheath_L, previous_sheath_final_hidden_L)

            if not prev_leaf_emerged: #: Before the emergence of the previous leaf. Exponential-like elong.
                ## delta leaf length
                delta_leaf_L = model.calculate_deltaL_preE(hiddenzone_inputs['sucrose'], hiddenzone_inputs['leaf_L'], hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'], self.delta_t)
                leaf_L = np.nanmin([curr_hiddenzone_outputs['leaf_Lmax'], (hiddenzone_inputs['leaf_L'] + delta_leaf_L)])

            else: #: After the emergence of the previous leaf.
                ## delta leaf length
                delta_leaf_L = model.calculate_deltaL_postE(hiddenzone_inputs['leaf_L'], curr_hiddenzone_outputs['leaf_Lmax'], hiddenzone_inputs['sucrose'], self.delta_t)
                leaf_L = np.nanmin([curr_hiddenzone_outputs['leaf_Lmax'], (hiddenzone_inputs['leaf_L'] + delta_leaf_L)])

                lamina_id = hiddenzone_id + tuple(['blade'])
                #: Lamina has not emerged
                if not curr_hiddenzone_outputs['leaf_is_emerged']:
                    #: Test of leaf emergence against hidden zone length. Assumes that a leaf cannot emerge before the previous one # TODO: besoin correction pour savoir à quel pas de temps exact??
                    curr_hiddenzone_outputs['leaf_is_emerged'] = model.calculate_leaf_emergence(hiddenzone_inputs['leaf_L'], curr_hiddenzone_outputs['hiddenzone_L'])
                    if curr_hiddenzone_outputs['leaf_is_emerged']: # Initialise lamina outputs
                        all_organs_outputs[lamina_id] = dict.fromkeys(ORGAN_OUTPUTS, 0)
                        all_organs_outputs[lamina_id]['is_growing'] = True

                        # Initialise variables for the next hidden zone
                        next_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] + 1])
                        if next_hiddenzone_id in all_hiddenzone_inputs:
                            next_hiddenzone_inputs = all_hiddenzone_inputs[next_hiddenzone_id]
                            next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]
                            next_hiddenzone_outputs['leaf_Lem_prev'] = next_hiddenzone_inputs['leaf_L']                                                                  # Leaf length at the time of the emergence of the previous leaf
                            next_hiddenzone_outputs['leaf_Lmax'] = model.calculate_leaf_Lmax(next_hiddenzone_outputs['leaf_Lem_prev'])                                   # Final leaf length
                            sheath_lamina_ratio = model.calculate_SL_ratio(next_hiddenzone_id[2])                                                                 # Sheath:Lamina final length ratio
                            next_hiddenzone_outputs['lamina_Lmax'] = model.calculate_lamina_Lmax(next_hiddenzone_outputs['leaf_Lmax'], sheath_lamina_ratio)              # Final lamina length
                            next_hiddenzone_outputs['sheath_Lmax'] = model.calculate_sheath_Lmax(next_hiddenzone_outputs['leaf_Lmax'], next_hiddenzone_outputs['lamina_Lmax'])  # Final sheath length
                            next_hiddenzone_outputs['leaf_Wmax'] = model.calculate_leaf_Wmax(next_hiddenzone_outputs['lamina_Lmax'], next_hiddenzone_inputs['fructan'], next_hiddenzone_inputs['mstruct'])        # Maximal leaf width
                            next_hiddenzone_outputs['SSLW'] = model.calculate_SSLW(next_hiddenzone_inputs['fructan'], next_hiddenzone_inputs['mstruct'] )                   # Structural Specific Lamina Weight
                            next_hiddenzone_outputs['SSSW'] = model.calculate_SSSW(next_hiddenzone_outputs['SSLW'])                                                      # Structural Specific Sheath Weight
                            self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                        else:
                            warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                #: Lamina has emerged and is growing
                elif curr_hiddenzone_outputs['leaf_is_emerged'] and all_organs_inputs[lamina_id]['is_growing']:
                    curr_organ_outputs = all_organs_outputs[lamina_id]
                    ## Length of emerged lamina
                    lamina_L = min(model.calculate_lamina_L(leaf_L, curr_hiddenzone_outputs['hiddenzone_L']), curr_hiddenzone_outputs['lamina_Lmax'])
                    curr_organ_outputs['visible_length'] = lamina_L
                    curr_organ_outputs['length'] = lamina_L

                    # Test end of elongation
                    if lamina_L >= curr_hiddenzone_outputs['lamina_Lmax']:
                        curr_organ_outputs['is_growing'] = False
                        curr_organ_outputs['final_hidden_length'] = 0
                        curr_organ_outputs['visible_length'] = lamina_L
                        curr_organ_outputs['length'] = lamina_L
                        # Initialise sheath outputs
                        sheath_id = hiddenzone_id + tuple(['sheath'])
                        all_organs_outputs[sheath_id] = dict.fromkeys(ORGAN_OUTPUTS, 0)
                        all_organs_outputs[sheath_id]['is_growing'] = True

                    # Update of lamina outputs
                    self.outputs['organs'][lamina_id] = curr_organ_outputs

                # Mature lamina, growing sheath
                else:
                    sheath_id = hiddenzone_id + tuple(['sheath'])
                    curr_organ_outputs = all_organs_outputs[sheath_id]

                    ## Length of emerged sheath
                    lamina_L = self.outputs['organs'][hiddenzone_id + tuple(['blade'])]['visible_length']
                    sheath_L = min(model.calculate_sheath_L(leaf_L, curr_hiddenzone_outputs['hiddenzone_L'], lamina_L), curr_hiddenzone_outputs['sheath_Lmax'])
                    curr_organ_outputs['visible_length'] = sheath_L

                    #: Test end of elong
                    if hiddenzone_inputs['leaf_L'] >= curr_hiddenzone_outputs['leaf_Lmax']: #TODO:  hiddenzone_inputs['leaf_L'] ou  hiddenzone_outputs['leaf_L']
                        curr_organ_outputs['final_hidden_length'] = curr_hiddenzone_outputs['hiddenzone_L']
                        curr_organ_outputs['length'] = sheath_L + curr_organ_outputs['final_hidden_length'] #: Length of the mature sheath = visible length + hiddenzone length
                        curr_organ_outputs['is_growing'] = False
                        curr_hiddenzone_outputs['leaf_is_growing'] = False

                    # Update of sheath outputs
                    self.outputs['organs'][sheath_id] = curr_organ_outputs

            # Update of leaf outputs, TODO: attention aux valeurs negatives
            curr_hiddenzone_outputs['leaf_L'] = leaf_L
            curr_hiddenzone_outputs['delta_leaf_L'] = np.nanmin([delta_leaf_L, (curr_hiddenzone_outputs['leaf_Lmax'] - hiddenzone_inputs['leaf_L'])])

            if curr_hiddenzone_outputs['leaf_is_growing']:
                self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs
            else: # End of leaf elong
                del self.outputs['hiddenzone'][hiddenzone_id]