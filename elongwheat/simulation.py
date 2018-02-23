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

import warnings
import copy

import numpy as np
import pandas as pd
import math

import model, parameters



#: the inputs needed by ElongWheat
HIDDENZONE_INPUTS = ['leaf_is_growing', 'internode_is_growing','leaf_pseudo_age','internode_pseudo_age','leaf_pseudostem_length','internode_distance_to_emergence', 'leaf_L', 'internode_L','leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'internode_Lmax','SSINW','internode_is_visible','sucrose', 'amino_acids', 'fructan','proteins', 'leaf_enclosed_mstruct','leaf_enclosed_Nstruct','internode_mstruct','internode_Nstruct','mstruct']
ELEMENT_INPUTS = ['length', 'is_growing']
SAM_INPUTS = ['sum_TT', 'status', 'nb_leaves', 'GA', 'height']

#: the outputs computed by ElongWheat # TODO : add be default all the attributes of the class HiddenZoneInit and ElementInit, and define which attribute is set by growthwheat.parameters or elongwheat.parameters
HIDDENZONE_OUTPUTS = ['leaf_is_growing', 'internode_is_growing','leaf_pseudo_age','internode_pseudo_age','leaf_pseudostem_length', 'delta_leaf_pseudostem_length','internode_distance_to_emergence', 'delta_internode_distance_to_emergence','leaf_L', 'delta_leaf_L', 'internode_L','delta_internode_L','leaf_Lmax', 'lamina_Lmax', 'sheath_Lmax', 'leaf_Wmax', 'SSLW', 'SSSW', 'leaf_is_emerged', 'internode_Lmax','SSINW','internode_is_visible', 'sucrose', 'amino_acids', 'fructan', 'proteins','leaf_enclosed_mstruct','leaf_enclosed_Nstruct','internode_mstruct','internode_Nstruct','mstruct']
ELEMENT_OUTPUTS = ['length','is_growing', 'diameter', 'sucrose', 'amino_acids','fructan','proteins','mstruct','Nstruct']
SAM_OUTPUTS = ['sum_TT', 'status', 'nb_leaves', 'GA', 'height']

#: the inputs and outputs of ElongWheat.
HIDDENZONE_INPUTS_OUTPUTS = sorted(set(HIDDENZONE_INPUTS + HIDDENZONE_OUTPUTS))
ELEMENT_INPUTS_OUTPUTS = sorted(set(ELEMENT_INPUTS + ELEMENT_OUTPUTS))
SAM_INPUTS_OUTPUTS = sorted(set(SAM_INPUTS + SAM_OUTPUTS))

#: topology colums for ligule height dataframe
LIGULE_TOPOLOGY_COLUMNS = ['SAM_id', 'phytomer', 'ligule height']

class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass


class Simulation(object):
    """The Simulation class allows to initialize and run a simulation.
    """

    def __init__(self, delta_t=1):

        #: The inputs of elong-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element): {element_input_name: element_input_value, ...}, ...},
        #:      'SAM': {(plant_index, axis_label): {SAM_input_name: SAM_input_value, ...}, ...},
        #:      'sheath_internode_lengths': {(plant_index, axis_label, metamer_index): {'sheath': [list of sheath length belonging to the phytomer],
        #:                                                                             'cumulated_internode': [list of internode lengths cumulated from phytomer 1 to n]}, ...}
        #:
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of elong-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element): {element_output_name: element_output_value, ...}, ...},
        #:      'SAM': {(plant_index, axis_label): {SAM_output_name: SAM_output_value, ...}, ...}}
        #:
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

    def run(self, Tair, Tsol):
        """
        Run the simulation.

        :Parameters:

            - `Tair` (:class:`float`) - air temperature at t (degree Celsius)
            - `Tsol` (:class:`float`) - soil temperature at t (degree Celsius)
        """

        # Copy the inputs into the output dict
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hiddenzone', 'elements', 'SAM', 'sheath_internode_lengths'])})

        # Hidden zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # Elements
        all_element_inputs = self.inputs['elements']
        all_element_outputs = self.outputs['elements']

        # SAM
        all_SAM_inputs = self.inputs['SAM']
        all_SAM_outputs = self.outputs['SAM']

        # sheath and internode lengths
        all_sheath_internode_lengths = self.inputs['sheath_internode_lengths']

        # Ligule heights
        all_ligule_height_df = pd.DataFrame(columns=LIGULE_TOPOLOGY_COLUMNS)

        # --  Beginning of computations

        # SAM
        for SAM_id, SAM_inputs in sorted(all_SAM_inputs.iteritems()):
            curr_SAM_outputs = all_SAM_outputs[SAM_id]
            nb_leaves = curr_SAM_outputs['nb_leaves']

            ## height of the SAM
            below_internode_lengths = all_sheath_internode_lengths[SAM_id][nb_leaves]['cumulated_internode']
            SAM_height = model.calculate_cumulated_internode_length(below_internode_lengths)
            curr_SAM_outputs['height'] = SAM_height

            ## temperature choice
            growth_temperature = model.calculate_growing_temperature(Tair, Tsol, SAM_height)

            ## update sum_TT
            curr_SAM_outputs['sum_TT'] = model.calculate_SAM_sumTT(growth_temperature, SAM_inputs['sum_TT'], self.delta_t)

            ## update SAM status, leaf number and
            init_leaf, curr_SAM_outputs['nb_leaves'], curr_SAM_outputs['status'] = model.calculate_SAM_status(SAM_inputs['status'], curr_SAM_outputs['sum_TT'], nb_leaves)

            ## GA production
            curr_SAM_outputs['GA'] = model.calculate_SAM_GA(curr_SAM_outputs['status'], curr_SAM_outputs['sum_TT'])

            ## add hiddenzone
            for i in range(0, init_leaf):
                # Initialise hiddenzone
                hiddenzone_id = SAM_id + tuple([1 + i + curr_SAM_outputs['nb_leaves'] - init_leaf]) # TODO: peut etre simplifié tant que 'calculate_SAM_status' renvoie 1 erreur dans le cas init_leaf>1
                new_hiddenzone = parameters.HiddenZoneInit().__dict__
                self.outputs['hiddenzone'][hiddenzone_id] = new_hiddenzone

            ## Ligule height
            all_ligule_height_df = model.calculate_ligule_height(all_sheath_internode_lengths[SAM_id], all_element_inputs, SAM_id, all_ligule_height_df)

            self.outputs['SAM'][SAM_id] = curr_SAM_outputs

        # Hiddenzones
        for hiddenzone_id, hiddenzone_inputs in sorted(all_hiddenzone_inputs.iteritems()):
            SAM_id = hiddenzone_id[:2]
            phytomer_id = hiddenzone_id[2]
            curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]
            curr_SAM_outputs = all_SAM_outputs[SAM_id]

            # Found previous hidden zone
            prev_hiddenzone_id = tuple(list(SAM_id) + [phytomer_id - 1])
            if prev_hiddenzone_id in all_hiddenzone_inputs:
                prev_leaf_emerged = all_hiddenzone_inputs[prev_hiddenzone_id]['leaf_is_emerged']
            else:
                prev_leaf_emerged = True

            # Cumulated length of internodes up to the hidden zone
            below_internode_lengths = all_sheath_internode_lengths[SAM_id][phytomer_id]['cumulated_internode']
            bottom_hiddenzone_height = model.calculate_cumulated_internode_length(below_internode_lengths)

            # Distance between the bottom of the hiddenzone and the highest previous ligule
            leaf_pseudostem_length = model.calculate_leaf_pseudostem_length(all_ligule_height_df[all_ligule_height_df['SAM_id']==SAM_id], bottom_hiddenzone_height, phytomer_id)
            curr_hiddenzone_outputs['leaf_pseudostem_length'] =  leaf_pseudostem_length
            curr_hiddenzone_outputs['delta_leaf_pseudostem_length'] = leaf_pseudostem_length - hiddenzone_inputs['leaf_pseudostem_length'] # Cette variable est utilisee dans growthwheat

            # Calculate the internode pseudostem length
            curr_internode_L = hiddenzone_inputs['internode_L']
            internode_distance_to_emergence = model.calculate_internode_distance_to_emergence(leaf_pseudostem_length, curr_internode_L)
            curr_hiddenzone_outputs['internode_distance_to_emergence'] = internode_distance_to_emergence
            curr_hiddenzone_outputs['delta_internode_distance_to_emergence'] = internode_distance_to_emergence - hiddenzone_inputs['internode_distance_to_emergence']  # Cette variable est utilisee dans growthwheat

            # In case leaf is already mature but internode is growing, we update sheath visible and hidden lengths.
            if not curr_hiddenzone_outputs['leaf_is_growing'] and curr_hiddenzone_outputs['leaf_is_emerged']:
               hidden_sheath_id = hiddenzone_id + tuple(['sheath', 'HiddenElement'])
               visible_sheath_id = hiddenzone_id + tuple(['sheath', 'StemElement'])
               total_sheath_L = self.inputs['elements'][hidden_sheath_id]['length'] + self.inputs['elements'][visible_sheath_id]['length']
               updated_sheath_hidden_length = min(total_sheath_L, leaf_pseudostem_length)
               updated_sheath_visible_length = total_sheath_L - updated_sheath_hidden_length
               self.outputs['elements'][hidden_sheath_id]['length'] = updated_sheath_hidden_length
               self.outputs['elements'][visible_sheath_id]['length'] = updated_sheath_visible_length

            #: Leaf elongation
            if curr_hiddenzone_outputs['leaf_is_growing']:

               if leaf_pseudostem_length < 0:
                  print 'Pseudostem length of {} decreased while leaf growing.'.format(hiddenzone_id)

               if not prev_leaf_emerged: #: Before the emergence of the previous leaf. Exponential-like elongation.
                  ## delta leaf length
                  delta_leaf_L = model.calculate_deltaL_preE(hiddenzone_inputs['sucrose'], hiddenzone_inputs['leaf_L'], hiddenzone_inputs['amino_acids'], hiddenzone_inputs['leaf_enclosed_mstruct'], self.delta_t, phytomer_id)
                  leaf_L = hiddenzone_inputs['leaf_L'] + delta_leaf_L

               else: #: After the emergence of the previous leaf.
                     ## delta leaf length
                    leaf_pseudo_age = model.calculate_leaf_pseudo_age(hiddenzone_inputs['leaf_pseudo_age'], self.delta_t)
                    curr_hiddenzone_outputs['leaf_pseudo_age'] = leaf_pseudo_age
                    leaf_L = model.calculate_L_postE(leaf_pseudo_age, curr_hiddenzone_outputs['leaf_Lmax'])
                    delta_leaf_L = leaf_L - hiddenzone_inputs['leaf_L']

                    lamina_id = hiddenzone_id + tuple(['blade', 'LeafElement1'])
                    #: Lamina has not emerged
                    if not curr_hiddenzone_outputs['leaf_is_emerged']:
                        #: Test of leaf emergence against distance to leaf emergence. Assumes that a leaf cannot emerge before the previous one # TODO: besoin correction pour savoir a quel pas de temps exact??
                        curr_hiddenzone_outputs['leaf_is_emerged'] = model.calculate_leaf_emergence(hiddenzone_inputs['leaf_L'], leaf_pseudostem_length)
                        if curr_hiddenzone_outputs['leaf_is_emerged']: # Initialise lamina outputs
                            new_lamina = parameters.ElementInit().__dict__
                            self.outputs['elements'][lamina_id] = new_lamina
                            ## Length of emerged lamina
                            curr_lamina_outputs = all_element_outputs[lamina_id]
                            lamina_L = min(model.calculate_lamina_L(leaf_L, leaf_pseudostem_length, hiddenzone_id), curr_hiddenzone_outputs['lamina_Lmax'])
                            curr_lamina_outputs['length'] = lamina_L
                            # Update of lamina outputs
                            self.outputs['elements'][lamina_id] = curr_lamina_outputs

                            # Initialise variables for the next hidden zone as its previous leaf has now emerged
                            next_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2] + 1])
                            if next_hiddenzone_id in all_hiddenzone_inputs:
                                next_hiddenzone_inputs = all_hiddenzone_inputs[next_hiddenzone_id]
                                next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]
                                next_hiddenzone_outputs['leaf_Lmax'] = model.calculate_leaf_Lmax(next_hiddenzone_inputs['leaf_L'])                                                  # Final leaf length
                                sheath_lamina_ratio = model.calculate_SL_ratio(next_hiddenzone_id[2])                                                                               # Sheath:Lamina final length ratio
                                next_hiddenzone_outputs['lamina_Lmax'] = model.calculate_lamina_Lmax(next_hiddenzone_outputs['leaf_Lmax'], sheath_lamina_ratio)                     # Final lamina length
                                next_hiddenzone_outputs['sheath_Lmax'] = model.calculate_sheath_Lmax(next_hiddenzone_outputs['leaf_Lmax'], next_hiddenzone_outputs['lamina_Lmax'])  # Final sheath length
                                next_hiddenzone_outputs['leaf_Wmax'] = model.calculate_leaf_Wmax(next_hiddenzone_outputs['lamina_Lmax'], next_hiddenzone_inputs['fructan'], next_hiddenzone_inputs['mstruct']) # Maximal leaf width
                                next_hiddenzone_outputs['SSLW'] = model.calculate_SSLW(next_hiddenzone_inputs['fructan'], next_hiddenzone_inputs['mstruct'] )                       # Structural Specific Lamina Weight
                                next_hiddenzone_outputs['SSSW'] = model.calculate_SSSW(next_hiddenzone_outputs['SSLW'])                                                             # Structural Specific Sheath Weight
                                next_hiddenzone_outputs['leaf_pseudo_age'] = 0                                                                                                      # Pseudo age of the leaf since beginning of automate growth (s)
                                self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                            else:
                                warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                    #: Lamina has emerged and is growing
                    elif curr_hiddenzone_outputs['leaf_is_emerged'] and all_element_inputs[lamina_id]['is_growing']:
                        curr_lamina_outputs = all_element_outputs[lamina_id]
                        ## Length of emerged lamina
                        lamina_L = min(model.calculate_lamina_L(leaf_L, leaf_pseudostem_length, hiddenzone_id), curr_hiddenzone_outputs['lamina_Lmax'])
                        curr_lamina_outputs['length'] = lamina_L

                        # Test end of elongation
                        if lamina_L >= curr_hiddenzone_outputs['lamina_Lmax']:
                            curr_lamina_outputs['is_growing'] = False
                            curr_lamina_outputs['length'] = curr_hiddenzone_outputs['lamina_Lmax']
                            # Initialise sheath outputs
                            visible_sheath_id = hiddenzone_id + tuple(['sheath', 'StemElement'])
                            new_sheath = parameters.ElementInit().__dict__
                            self.outputs['elements'][visible_sheath_id] = new_sheath
                            ## Length of emerged sheath
                            curr_visible_sheath_outputs = all_element_outputs[visible_sheath_id]
                            emerged_sheath_L = min(model.calculate_emerged_sheath_L(leaf_L, leaf_pseudostem_length, lamina_L), curr_hiddenzone_outputs['sheath_Lmax'])
                            curr_visible_sheath_outputs['length'] = emerged_sheath_L
                            # Update of sheath outputs
                            self.outputs['elements'][visible_sheath_id] = curr_visible_sheath_outputs

                            # Initialise variables for the next internode
                            next_hiddenzone_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2]+1])
                            if next_hiddenzone_id in all_hiddenzone_inputs:
                               next_hiddenzone_outputs = all_hiddenzone_outputs[next_hiddenzone_id]
                               if curr_SAM_outputs['GA']:
                                  next_hiddenzone_outputs['internode_Lmax'] = model.calculate_internode_Lmax(next_hiddenzone_outputs['internode_L'])      #: Final internode length
                               next_hiddenzone_outputs['SSIW'] = model.calculate_SSIW(next_hiddenzone_outputs['SSSW'])                                    #: Structural Specific Internode Weight
                               next_hiddenzone_outputs['internode_pseudo_age'] = 0                                                                        #: Pseudo age of the internode since beginning of automate growth (s)
                               self.outputs['hiddenzone'][next_hiddenzone_id] = next_hiddenzone_outputs
                            else:
                                warnings.warn('No next hidden zone found for hiddenzone {}.'.format(hiddenzone_id))

                        # Update of lamina outputs
                        self.outputs['elements'][lamina_id] = curr_lamina_outputs

                    # Mature lamina, growing sheath
                    else:
                        visible_sheath_id = hiddenzone_id + tuple(['sheath', 'StemElement'])
                        curr_visible_sheath_outputs = all_element_outputs[visible_sheath_id]

                        ## Length of emerged sheath
                        lamina_L = self.outputs['elements'][hiddenzone_id + tuple(['blade', 'LeafElement1'])]['length']
                        emerged_sheath_L = min(model.calculate_emerged_sheath_L(leaf_L, leaf_pseudostem_length, lamina_L), curr_hiddenzone_outputs['sheath_Lmax'])
                        curr_visible_sheath_outputs['length'] = emerged_sheath_L

                        #: Test end of elongation
                        if leaf_L >= curr_hiddenzone_outputs['leaf_Lmax']:
                            curr_visible_sheath_outputs['is_growing'] = False
                            curr_hiddenzone_outputs['leaf_is_growing'] = False
                            # Hidden sheath
                            hidden_sheath_id = hiddenzone_id + tuple(['sheath', 'HiddenElement'])
                            curr_hidden_sheath_outputs = dict.fromkeys(ELEMENT_OUTPUTS)
                            curr_hidden_sheath_outputs['length'] = leaf_pseudostem_length
                            curr_hidden_sheath_outputs['is_growing'] = False
                            self.outputs['elements'][hidden_sheath_id] = curr_hidden_sheath_outputs

                        # Update of sheath outputs
                        self.outputs['elements'][visible_sheath_id] = curr_visible_sheath_outputs


            #: Leaf is mature but internode may be growing
            else:
                 leaf_L = hiddenzone_inputs['leaf_L']
                 delta_leaf_L = 0

            # Update of leaf outputs
            curr_hiddenzone_outputs['leaf_L'] = leaf_L
            curr_hiddenzone_outputs['delta_leaf_L'] = delta_leaf_L

            #: Internode elongation
            ## Initialisation of internode elongation
            if (not curr_hiddenzone_outputs['internode_is_growing']) and curr_hiddenzone_outputs['internode_L'] == 0 :
               curr_hiddenzone_outputs['internode_is_growing'], curr_hiddenzone_outputs['internode_L'] = model.calculate_init_internode_elongation(curr_SAM_outputs['sum_TT'], phytomer_id) # As for leaf primordia, we neglect CN growth due to IN length initialisation

            if curr_hiddenzone_outputs['internode_is_growing']:

               # Found previous lamina to know if the previous leaf is ligulated.
               prev_lamina_id = tuple(list(hiddenzone_id[:2]) + [hiddenzone_id[2]-1]) + tuple(['blade', 'LeafElement1'])
               if prev_lamina_id in all_element_inputs:
                  prev_leaf_ligulated = not all_element_inputs[prev_lamina_id]['is_growing']
               else:
                  prev_leaf_ligulated = False

               #: Before ligulation of the leaf on the previous phytomer. Exponential-like elong. cf. Kirby 1988, Malvoisin 1984 II
               if not prev_leaf_ligulated:
                 delta_internode_L = model.calculate_delta_internode_L_preL(phytomer_id, curr_hiddenzone_outputs['sucrose'], curr_hiddenzone_outputs['internode_L'], curr_hiddenzone_outputs['amino_acids'], curr_hiddenzone_outputs['mstruct'], self.delta_t)
                 internode_L = curr_hiddenzone_outputs['internode_L'] + delta_internode_L  # TODO: Ckeck internode_L is not too large (in the case of long delta_t)

               #: After ligulation of the leaf on the previous phytomer.
               else:
                   internode_pseudo_age = model.calculate_internode_pseudo_age(curr_hiddenzone_outputs['internode_pseudo_age'], self.delta_t)
                   curr_hiddenzone_outputs['internode_pseudo_age'] = internode_pseudo_age

                   ## Elongation only if Gibberelin production by SAM
                   if curr_SAM_outputs['GA'] :

                       # Case of internodes that will not fully elongate, GA synthesis started after their previous leaf ligulation (ie no Lmax defined)
                       if np.isnan(curr_hiddenzone_outputs['internode_Lmax']):
                          curr_hiddenzone_outputs['internode_Lmax'] = model.calculate_short_internode_Lmax(curr_hiddenzone_outputs['internode_L'], curr_hiddenzone_outputs['internode_pseudo_age'])

                       internode_L = model.calculate_internode_L_postL(curr_hiddenzone_outputs['internode_pseudo_age'], curr_hiddenzone_outputs['internode_Lmax'])
                       delta_internode_L = internode_L - curr_hiddenzone_outputs['internode_L']

                       internode_id = hiddenzone_id + tuple(['internode', 'StemElement'])
                       #: Internode is not visible
                       if not curr_hiddenzone_outputs['internode_is_visible']:
                          #: Test of internode emergence.
                          curr_hiddenzone_outputs['internode_is_visible'] = model.calculate_internode_visibility(curr_hiddenzone_outputs['internode_L'], internode_distance_to_emergence)
                          if curr_hiddenzone_outputs['internode_is_visible']: # Initialise internode outputs
                             new_internode_outputs = parameters.ElementInit().__dict__
                             self.outputs['elements'][internode_id] = new_internode_outputs
                             self.outputs['elements'][internode_id]['length'] = min(curr_hiddenzone_outputs['internode_Lmax'], model.calculate_emerged_internode_L(internode_L, internode_distance_to_emergence))
                       #: Internode is visible
                       else:
                            self.outputs['elements'][internode_id]['length'] = min(curr_hiddenzone_outputs['internode_Lmax'], model.calculate_emerged_internode_L(internode_L, internode_distance_to_emergence))

                       # Test end of elongation
                       if model.calculate_end_internode_elongation(internode_L, curr_hiddenzone_outputs['internode_Lmax'], curr_hiddenzone_outputs['internode_pseudo_age']):
                          curr_hiddenzone_outputs['internode_is_growing'] = False
                          # Visible internode
                          if curr_hiddenzone_outputs['internode_is_visible']:
                             self.outputs['elements'][internode_id]['is_growing'] = False
                          # Hidden internode
                          hidden_internode_id = hiddenzone_id + tuple(['internode', 'HiddenElement'])
                          curr_hidden_internode_outputs = dict.fromkeys(ELEMENT_OUTPUTS)
                          curr_hidden_internode_outputs['length'] = min(internode_distance_to_emergence, curr_hiddenzone_outputs['internode_Lmax'])
                          curr_hidden_internode_outputs['is_growing'] = False
                          self.outputs['elements'][hidden_internode_id] = curr_hidden_internode_outputs

                   #: Internode not elongating because no GA
                   else:
                        internode_L = curr_hiddenzone_outputs['internode_L']
                        delta_internode_L = 0
                        # Test end of elongation
                        if model.calculate_end_internode_elongation(internode_L, curr_hiddenzone_outputs['internode_Lmax'], curr_hiddenzone_outputs['internode_pseudo_age']):
                          curr_hiddenzone_outputs['internode_is_growing'] = False
                          curr_hiddenzone_outputs['internode_Lmax'] = curr_hiddenzone_outputs['internode_L']
                          # Visible internode
                          if curr_hiddenzone_outputs['internode_is_visible']:
                             self.outputs['elements'][internode_id]['is_growing'] = False
                          # Hidden internode
                          hidden_internode_id = hiddenzone_id + tuple(['internode', 'HiddenElement'])
                          curr_hidden_internode_outputs = dict.fromkeys(ELEMENT_OUTPUTS)
                          curr_hidden_internode_outputs['length'] = min(internode_distance_to_emergence, curr_hiddenzone_outputs['internode_Lmax'])
                          curr_hidden_internode_outputs['is_growing'] = False
                          self.outputs['elements'][hidden_internode_id] = curr_hidden_internode_outputs

            #: Internode not elongating (not yet or already mature)
            else:
               internode_L = curr_hiddenzone_outputs['internode_L']
               delta_internode_L = 0

            # Update internodes outputs
            curr_hiddenzone_outputs['internode_L'] = internode_L
            curr_hiddenzone_outputs['delta_internode_L'] = delta_internode_L

            # Only growing hiddenzones are sent
            ## after end of elongation (leaf and/or internode), it should :
            ##       - pass by growth wheat for remobilisation
            ##       - pass another time by elong wheat for update of curr_element_outputs['final_hidden_length']
            ## the hiddenzone will then be deleted since both growing flags are False and both delta_L are zeros.
            ## For "internode_is_growing", the test is made on the inputs so we make sure it goes at least once for remobilisation by growthwheat.
            if hiddenzone_inputs['internode_is_growing'] or curr_hiddenzone_outputs['leaf_is_growing'] or curr_hiddenzone_outputs['delta_internode_L'] > 0 or curr_hiddenzone_outputs['delta_leaf_L'] > 0:
                self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs
            else: # End of internode elong
                del self.outputs['hiddenzone'][hiddenzone_id]