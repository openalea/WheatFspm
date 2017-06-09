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

import model, parameters
import warnings
import copy

import math ##MG

#: the inputs needed by ElongWheat
HIDDENZONE_INPUTS = ['leaf_is_growing', 'internode_is_growing','leaf_dist_to_emerge','internode_dist_to_emerge', 'leaf_L', 'internode_L','leaf_Lmax',  'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan','proteins', 'leaf_enclosed_mstruct','leaf_enclosed_Nstruct','internode_mstruct','internode_Nstruct','mstruct']
ORGAN_INPUTS = ['visible_length', 'is_growing', 'final_hidden_length', 'length']
SAM_INPUTS = ['sum_TT_prev_init', 'status','nb_leaves']

#: the outputs computed by ElongWheat
HIDDENZONE_OUTPUTS = ['leaf_is_growing', 'internode_is_growing','leaf_dist_to_emerge', 'delta_leaf_dist_to_emerge','internode_dist_to_emerge', 'delta_internode_dist_to_emerge','leaf_L', 'delta_leaf_L', 'internode_L','delta_internode_L','leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'sucrose', 'amino_acids', 'fructan', 'proteins','leaf_enclosed_mstruct','leaf_enclosed_Nstruct','internode_mstruct','internode_Nstruct','mstruct']
ORGAN_OUTPUTS = ['visible_length', 'is_growing', 'final_hidden_length', 'length']
SAM_OUTPUTS = ['sum_TT_prev_init', 'status','nb_leaves']

#: the inputs and outputs of ElongWheat.
HIDDENZONE_INPUTS_OUTPUTS = sorted(set(HIDDENZONE_INPUTS + HIDDENZONE_OUTPUTS))
ORGAN_INPUTS_OUTPUTS = sorted(set(ORGAN_INPUTS + ORGAN_OUTPUTS))
SAM_INPUTS_OUTPUTS = sorted(set(SAM_INPUTS + SAM_OUTPUTS))


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
        #:                                                                       'previous_sheath_final_hidden_length': previous_sheath_final_hidden_length,
        #:                                                                       'internode_length' : internode_length }, ...}}
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

    def run(self,Ta,Ts):
        """
        Run the simulation.

        :Parameters:

            - `Ta` (:class:`float`) - air temperature at t (degree Celsius)
            - `Ts` (:class:`float`) - soil temperature at t (degree Celsius)
        """

        # Copy the inputs into the output dict
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hiddenzone', 'organs', 'SAM','hiddenzone_L_calculation'])})

        # Hidden zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # organs
        all_organs_inputs = self.inputs['organs']
        all_organs_outputs = self.outputs['organs']

        # SAM
        all_SAM_inputs = self.inputs['SAM']
        all_SAM_outputs = self.outputs['SAM']

        # Temperature growing zone (used in SAM and hiddenzone)
        ## plant height
        plant_height=0
        for item in all_hiddenzone_inputs.values(): ##TODO: Find a nicer way to do it.
            plant_height+=item['internode_L']
        for organ_id, organ_inputs in all_organs_inputs.iteritems(): ##TODO: Find a nicer way to do it.
            if organ_id[3] == 'internode' and not organ_inputs['is_growing']:
               plant_height+=organ_inputs['length']

        ## temperature choice
        growth_temp = model.calculate_growing_temp(Ta, Ts, plant_height)


        # SAM
        for SAM_id, SAM_inputs in all_SAM_inputs.iteritems():
            curr_SAM_outputs = all_SAM_outputs[SAM_id]
            ## update sum_TT
            curr_SAM_outputs['sum_TT_prev_init'] , init_leaf = model.calculate_SAM_sumTT(growth_temp, SAM_inputs['sum_TT_prev_init'],  SAM_inputs['status'], self.delta_t)
            ## update SAM status
            curr_SAM_outputs['status'], curr_SAM_outputs['nb_leaves'] = model.calculate_SAM_status(SAM_inputs['status'],SAM_inputs['nb_leaves'],init_leaf)
            ## add hiddenzone
            for i in range(0,init_leaf):
                # Initialise sheath outputs
                hiddenzone_id = SAM_id + tuple([ 1+i+curr_SAM_outputs['nb_leaves']-init_leaf ])
                new_hiddenzone = parameters.HiddenZoneInit().__dict__
                self.outputs['hiddenzone'][hiddenzone_id] = new_hiddenzone

            self.outputs['SAM'][SAM_id] = curr_SAM_outputs

        # Hidden zone lengths
        all_hiddenzone_L_calculation_inputs = self.inputs['hiddenzone_L_calculation']

        for hiddenzone_id, hiddenzone_inputs in all_hiddenzone_inputs.iteritems():
            curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]
            # Found previous hidden zone TODO: a améliorer
            prev_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] - 1])
            if prev_hiddenzone_id in all_hiddenzone_inputs:
                prev_leaf_emerged = all_hiddenzone_inputs[prev_hiddenzone_id]['leaf_is_emerged']
            else:
                prev_leaf_emerged = True

            # Get values for update distance for the leaf to emerge and distance for internode to be visible
            previous_hiddenzone_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['previous_leaf_dist_to_emerge']
            curr_internode_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['internode_length']
            previous_sheath_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['previous_sheath_visible_length']
            previous_sheath_final_hidden_L = all_hiddenzone_L_calculation_inputs[hiddenzone_id]['previous_sheath_final_hidden_length']

            # Update distance for the leaf to emerge
            curr_hiddenzone_outputs['leaf_dist_to_emerge'] = model.calculate_leaf_dist_to_emerge(previous_hiddenzone_L, curr_internode_L, previous_sheath_L, previous_sheath_final_hidden_L)
            curr_hiddenzone_outputs['delta_leaf_dist_to_emerge'] = curr_hiddenzone_outputs['leaf_dist_to_emerge'] - hiddenzone_inputs['leaf_dist_to_emerge']

            # In case leaf is already mature but internode is growing, we update sheath hidden length.
            if not curr_hiddenzone_outputs['leaf_is_growing'] and curr_hiddenzone_outputs['leaf_is_emerged']:
               sheath_id = hiddenzone_id + tuple(['sheath'])
               self.outputs['organs'][sheath_id]['final_hidden_length'] = curr_hiddenzone_outputs['leaf_dist_to_emerge']
               self.outputs['organs'][sheath_id]['visible_length'] = self.outputs['organs'][sheath_id]['length'] - curr_hiddenzone_outputs['leaf_dist_to_emerge']

            # Update distance for the internode to be visible
            curr_hiddenzone_outputs['internode_dist_to_emerge'] = model.calculate_internode_dist_to_emerge(previous_hiddenzone_L, previous_sheath_L, previous_sheath_final_hidden_L)
            curr_hiddenzone_outputs['delta_internode_dist_to_emerge'] = curr_hiddenzone_outputs['internode_dist_to_emerge'] - hiddenzone_inputs['internode_dist_to_emerge']


            #: Leaf elongation
            if curr_hiddenzone_outputs['leaf_is_growing']:

               if not prev_leaf_emerged: #: Before the emergence of the previous leaf. Exponential-like elong.
                  ## delta leaf length
                  delta_leaf_L = model.calculate_deltaL_preE(hiddenzone_inputs['sucrose'], hiddenzone_inputs['leaf_L'], hiddenzone_inputs['amino_acids'], hiddenzone_inputs['leaf_enclosed_mstruct'], self.delta_t)
                  leaf_L = hiddenzone_inputs['leaf_L'] + delta_leaf_L # TODO: Ckeck leaf_L is not too large (in the case of long delta_t) ##MG : it was np.nanmin([curr_hiddenzone_outputs['leaf_Lmax'], (hiddenzone_inputs['leaf_L'] + delta_leaf_L)]) but leaf_Lmax not defined yet

               else: #: After the emergence of the previous leaf.
                     ## delta leaf length
                    delta_leaf_L = model.calculate_deltaL_postE(hiddenzone_inputs['leaf_L'], curr_hiddenzone_outputs['leaf_Lmax'], hiddenzone_inputs['sucrose'], self.delta_t)
                    leaf_L = np.nanmin([curr_hiddenzone_outputs['leaf_Lmax'], (hiddenzone_inputs['leaf_L'] + delta_leaf_L)])

                    lamina_id = hiddenzone_id + tuple(['blade'])
                    #: Lamina has not emerged
                    if not curr_hiddenzone_outputs['leaf_is_emerged']:
                        #: Test of leaf emergence against distance to leaf emergence. Assumes that a leaf cannot emerge before the previous one # TODO: besoin correction pour savoir à quel pas de temps exact??
                        curr_hiddenzone_outputs['leaf_is_emerged'] = model.calculate_leaf_emergence(hiddenzone_inputs['leaf_L'], curr_hiddenzone_outputs['leaf_dist_to_emerge']) # leaf_L from previous step (hiddenzone_inputs['leaf_L']) should be tested against current calculation of pseudostem length (curr_hiddenzone_outputs['leaf_dist_to_emerge']) for constitency
                        if curr_hiddenzone_outputs['leaf_is_emerged']: # Initialise lamina outputs
                            new_lamina_outputs = dict.fromkeys(ORGAN_OUTPUTS, 0)
                            new_lamina_outputs['is_growing'] = True
                            self.outputs['organs'][lamina_id] = new_lamina_outputs

                            # Initialise variables for the next hidden zone
                            next_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] + 1])
                            if next_hiddenzone_id in all_hiddenzone_inputs:
                                next_hiddenzone_inputs = all_hiddenzone_inputs[next_hiddenzone_id]
                                next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]
                                next_hiddenzone_outputs['leaf_Lmax'] = model.calculate_leaf_Lmax(next_hiddenzone_inputs['leaf_L'])                                   # Final leaf length
                                sheath_lamina_ratio = model.calculate_SL_ratio(next_hiddenzone_id[2])                                                                 # Sheath:Lamina final length ratio
                                next_hiddenzone_outputs['lamina_Lmax'] = model.calculate_lamina_Lmax(next_hiddenzone_outputs['leaf_Lmax'], sheath_lamina_ratio)              # Final lamina length
                                next_hiddenzone_outputs['sheath_Lmax'] = model.calculate_sheath_Lmax(next_hiddenzone_outputs['leaf_Lmax'], next_hiddenzone_outputs['lamina_Lmax'])  # Final sheath length
                                next_hiddenzone_outputs['leaf_Wmax'] = model.calculate_leaf_Wmax(next_hiddenzone_outputs['lamina_Lmax'], next_hiddenzone_inputs['fructan'], next_hiddenzone_inputs['mstruct'])        # Maximal leaf width
                                next_hiddenzone_outputs['SSLW'] = model.calculate_SSLW(next_hiddenzone_inputs['fructan'], next_hiddenzone_inputs['leaf_enclosed_mstruct'] )                   # Structural Specific Lamina Weight
                                next_hiddenzone_outputs['SSSW'] = model.calculate_SSSW(next_hiddenzone_outputs['SSLW'])                                                      # Structural Specific Sheath Weight
                                self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                            else:
                                warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                    #: Lamina has emerged and is growing
                    elif curr_hiddenzone_outputs['leaf_is_emerged'] and all_organs_inputs[lamina_id]['is_growing']:
                        curr_organ_outputs = all_organs_outputs[lamina_id]
                        ## Length of emerged lamina
                        lamina_L = min(model.calculate_lamina_L(leaf_L, curr_hiddenzone_outputs['leaf_dist_to_emerge']), curr_hiddenzone_outputs['lamina_Lmax'])
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
                            new_sheath = dict.fromkeys(ORGAN_OUTPUTS, 0)
                            new_sheath[sheath_id]['is_growing'] = True
                            self.outputs['organs'][sheath_id] = new_sheath

                        # Update of lamina outputs
                        self.outputs['organs'][lamina_id] = curr_organ_outputs

                    # Mature lamina, growing sheath
                    else:
                        sheath_id = hiddenzone_id + tuple(['sheath'])
                        curr_organ_outputs = all_organs_outputs[sheath_id]

                        ## Length of emerged sheath
                        lamina_L = self.outputs['organs'][hiddenzone_id + tuple(['blade'])]['visible_length']
                        sheath_L = min(model.calculate_sheath_L(leaf_L, curr_hiddenzone_outputs['leaf_dist_to_emerge'], lamina_L), curr_hiddenzone_outputs['sheath_Lmax'])
                        curr_organ_outputs['visible_length'] = sheath_L

                        #: Test end of elong
                        if leaf_L >= curr_hiddenzone_outputs['leaf_Lmax']: ## Used to be hiddenzone_inputs['leaf_L'] instead of leaf_L
                            curr_organ_outputs['final_hidden_length'] = curr_hiddenzone_outputs['leaf_dist_to_emerge']
                            curr_organ_outputs['length'] = sheath_L + curr_organ_outputs['final_hidden_length'] #: Length of the mature sheath = visible length + pseudostem length
                            curr_organ_outputs['is_growing'] = False
                            curr_hiddenzone_outputs['leaf_is_growing'] = False
                            curr_hiddenzone_outputs['internode_is_growing'] = True

                        # Update of sheath outputs
                        self.outputs['organs'][sheath_id] = curr_organ_outputs

            #: Leaf is mature
            else:
                 leaf_L = hiddenzone_inputs['leaf_L']
                 delta_leaf_L = 0

            # Update of leaf outputs, #TODO: attention aux valeurs negatives
            curr_hiddenzone_outputs['leaf_L'] = leaf_L

            if curr_hiddenzone_outputs['leaf_Lmax'] == None or math.isnan(curr_hiddenzone_outputs['leaf_Lmax']): #TODO: Find a nicer way
               curr_hiddenzone_outputs['delta_leaf_L'] = delta_leaf_L
            else:
                curr_hiddenzone_outputs['delta_leaf_L'] = np.nanmin([delta_leaf_L, (curr_hiddenzone_outputs['leaf_Lmax'] - hiddenzone_inputs['leaf_L'])])

            #: Internode elongation
            if curr_hiddenzone_outputs['internode_is_growing']:

               delta_internode_L = model.calculate_delta_internode_L(hiddenzone_inputs['internode_L'], self.delta_t)
               internode_L = np.nanmin([curr_hiddenzone_outputs['internode_dist_to_emerge'], (hiddenzone_inputs['internode_L'] + delta_internode_L)])

               # Test end of elongation
               if model.calculate_end_internode_elongation(internode_L, curr_hiddenzone_outputs['internode_dist_to_emerge']): #TODO: add visible part of IN
                  curr_hiddenzone_outputs['internode_is_growing'] = False
                  # Initialise internode outputs
                  internode_id = hiddenzone_id + tuple(['internode'])
                  new_internode = dict.fromkeys(ORGAN_OUTPUTS, 0)
                  new_internode['is_growing'] = False
                  new_internode['final_hidden_length'] = internode_L #TODO: hiddenzone_inputs['internode_L'] ou internode_L
                  new_internode['length'] = internode_L #TODO: hiddenzone_inputs['internode_L'] ou internode_L
                  new_internode['visible_length'] = 0
                  self.outputs['organs'][internode_id] = new_internode

            #: IN not elongating (not yet or already mature)
            else:
               internode_L = hiddenzone_inputs['internode_L']
               delta_internode_L = 0


            # Update internodes outputs
            curr_hiddenzone_outputs['internode_L'] = internode_L
            curr_hiddenzone_outputs['delta_internode_L'] = np.nanmin([delta_internode_L, (curr_hiddenzone_outputs['internode_dist_to_emerge'] - hiddenzone_inputs['internode_L'])])

            # Only growing hiddenzones are sent
            ## after end of elongation (leaf and/or internode), it should :
            ##       - pass by growth wheat for CN flux + remobilisation
            ##       - pass another time by elong wheat for update of curr_organ_outputs['final_hidden_length']
            ## the hiddenzone will then be deleted since both growing flags are False and both delta_L are zeros.

            if (curr_hiddenzone_outputs['internode_is_growing'] or curr_hiddenzone_outputs['leaf_is_growing'] or curr_hiddenzone_outputs['delta_internode_L'] > 0 or curr_hiddenzone_outputs['delta_leaf_L'] > 0):
                self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs
            else: # End of internode elong
                del self.outputs['hiddenzone'][hiddenzone_id]
