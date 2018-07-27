# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    growthwheat.simulation
    ~~~~~~~~~~~~~~~~~~

    The module :mod:`growthwheat.simulation`.

    :copyright: Copyright 2014-2015 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2015.
"""

"""
    Information about this versioned file:
        $LastChangedBy: mngauthier $
        $LastChangedDate: 2018-07-20 15:10:01 +0200 (ven., 20 juil. 2018) $
        $LastChangedRevision: 40 $
        $URL: https://subversion.renater.fr/growth-wheat/trunk/trunk/growthwheat/simulation.py $
        $Id: simulation.py 40 2018-07-20 13:10:01Z mngauthier $
"""

import model, parameters
import copy

from respiwheat.model import RespirationModel

#: the inputs needed by GrowthWheat
HIDDENZONE_INPUTS = ['leaf_is_growing','internode_is_growing','leaf_L', 'delta_leaf_L', 'internode_L','delta_internode_L', 'leaf_pseudostem_length', 'delta_leaf_pseudostem_length',
                     'internode_distance_to_emergence','delta_internode_distance_to_emergence', 'SSLW', 'LSSW', 'LSIW','leaf_is_emerged', 'internode_is_visible','leaf_pseudo_age','internode_pseudo_age',
                     'sucrose', 'amino_acids', 'fructan','proteins','leaf_enclosed_mstruct','leaf_enclosed_Nstruct','internode_enclosed_mstruct', 'internode_enclosed_Nstruct', 'mstruct']
ELEMENT_INPUTS = ['is_growing', 'mstruct', 'green_area','length', 'sucrose', 'amino_acids', 'fructan','proteins','Nstruct']
ROOT_INPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct']

#: the outputs computed by GrowthWheat
HIDDENZONE_OUTPUTS = ['sucrose', 'amino_acids','fructan','proteins', 'leaf_enclosed_mstruct', 'leaf_enclosed_Nstruct', 'internode_enclosed_mstruct', 'internode_enclosed_Nstruct','mstruct', 'Respi_growth', 'sucrose_consumption_mstruct', 'AA_consumption_mstruct']
ELEMENT_OUTPUTS = ['sucrose', 'amino_acids','fructan','proteins', 'mstruct', 'Nstruct','green_area']
ROOT_OUTPUTS = ['sucrose', 'amino_acids', 'mstruct', 'Nstruct', 'Respi_growth', 'rate_mstruct_growth']

#: the inputs and outputs of GrowthWheat.
HIDDENZONE_INPUTS_OUTPUTS = sorted(set(HIDDENZONE_INPUTS + HIDDENZONE_OUTPUTS))
ELEMENT_INPUTS_OUTPUTS = sorted(set(ELEMENT_INPUTS + ELEMENT_OUTPUTS))
ROOT_INPUTS_OUTPUTS = sorted(set(ROOT_INPUTS + ROOT_OUTPUTS))


class SimulationError(Exception): pass
class SimulationRunError(SimulationError): pass


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self, delta_t=1):

        #: The inputs of growth-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {organ_input_name: organ_input_value, ...}, ...},
        #:      'roots': {(plant_index, axis_label): {root_input_name: root_input_value, ...}, ...}}
        #: See :TODO?
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of growth-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'hiddenzone': {(plant_index, axis_label, metamer_index): {hiddenzone_input_name: hiddenzone_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {organ_input_name: organ_input_value, ...}, ...}
        #:      'roots': {(plant_index, axis_label): {root_input_name: root_input_value, ...}, ...}}
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
        self.outputs.update({inputs_type: copy.deepcopy(all_inputs) for inputs_type, all_inputs in self.inputs.iteritems() if inputs_type in set(['hiddenzone', 'elements', 'roots'])})

        # Hidden growing zones
        all_hiddenzone_inputs = self.inputs['hiddenzone']
        all_hiddenzone_outputs = self.outputs['hiddenzone']

        # elements
        all_elements_inputs = self.inputs['elements']
        all_elements_outputs = self.outputs['elements']

        # roots
        all_roots_inputs = self.inputs['roots']
        all_roots_outputs = self.outputs['roots']

        # hidden zones and elements
        for hiddenzone_id, hiddenzone_inputs in all_hiddenzone_inputs.iteritems():

            curr_hiddenzone_outputs = all_hiddenzone_outputs[hiddenzone_id]

            # Initialisation of the exports towards the growing lamina or sheath
            delta_lamina_mstruct, delta_sheath_mstruct, delta_lamina_Nstruct, delta_sheath_Nstruct, export_sucrose, export_amino_acids = 0, 0, 0, 0, 0, 0

            # Delta Growth internode

            if not hiddenzone_inputs['internode_is_visible']: #: Internode is not visible
               ## delta mstruct of the internode
               delta_internode_enclosed_mstruct = model.calculate_delta_internode_enclosed_mstruct(hiddenzone_inputs['internode_L'], hiddenzone_inputs['delta_internode_L'])
               ## delta Nstruct of the internode
               delta_internode_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_internode_enclosed_mstruct)
            else : #: Internode is visible
                ## delta mstruct of the enclosed internode
                delta_internode_enclosed_mstruct = model.calculate_delta_internode_enclosed_mstruct(hiddenzone_inputs['internode_distance_to_emerge'], hiddenzone_inputs['delta_internode_distance_to_emerge'])
                ## delta Nstruct of the enclosed internode
                delta_internode_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_internode_enclosed_mstruct)

                visible_internode_id = hiddenzone_id + tuple(['internode','StemElement'])
                curr_visible_internode_inputs = all_elements_inputs[visible_internode_id]
                curr_visible_internode_outputs = all_elements_outputs[visible_internode_id]
                ## Delta mstruct of the emerged internode
                delta_internode_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['LSIW'], curr_visible_internode_inputs['mstruct'], curr_visible_internode_inputs['length'])
                ## Delta Nstruct of the emerged internode
                delta_internode_Nstruct = model.calculate_delta_Nstruct(delta_internode_mstruct)
                ## Export of sucrose from hiddenzone towards emerged internode
                export_sucrose = model.calculate_export(delta_internode_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                ## Export of amino acids from hiddenzone towards emerged internode
                export_amino_acids = model.calculate_export(delta_internode_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])

               # Update of internode outputs
                curr_visible_internode_outputs['mstruct'] += delta_internode_mstruct
                curr_visible_internode_outputs['Nstruct'] += delta_internode_Nstruct
                curr_visible_internode_outputs['sucrose'] += export_sucrose
                curr_visible_internode_outputs['amino_acids'] += export_amino_acids
                self.outputs['elements'][visible_internode_id] = curr_visible_internode_outputs

            # Delta Growth leaf

            if not hiddenzone_inputs['leaf_is_emerged']: #: leaf has not emerged
                ## delta mstruct of the hidden leaf
                delta_leaf_enclosed_mstruct = model.calculate_delta_leaf_enclosed_mstruct(hiddenzone_inputs['leaf_L'], hiddenzone_inputs['delta_leaf_L'])
                ## delta Nstruct of the hidden leaf
                delta_leaf_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_leaf_enclosed_mstruct)

            else: #: leaf has emerged
                ## delta mstruct of the enclosed leaf (which length is assumed to equal the length of the pseudostem)
                delta_leaf_enclosed_mstruct = model.calculate_delta_leaf_enclosed_mstruct_postE(hiddenzone_inputs['leaf_pseudostem_length'],
                                                                                                hiddenzone_inputs['delta_leaf_pseudostem_length'],
                                                                                                hiddenzone_inputs['mstruct'],
                                                                                                hiddenzone_inputs['LSSW'],
                                                                                                hiddenzone_inputs['leaf_pseudo_age'],
                                                                                                self.delta_t)
                ## delta Nstruct of the enclosed en leaf
                delta_leaf_enclosed_Nstruct = model.calculate_delta_Nstruct(delta_leaf_enclosed_mstruct)

                visible_lamina_id = hiddenzone_id + tuple(['blade','LeafElement1'])
                #: Lamina is growing
                if all_elements_inputs[visible_lamina_id]['is_growing']:
                    curr_visible_lamina_inputs = all_elements_inputs[visible_lamina_id]
                    curr_visible_lamina_outputs = all_elements_outputs[visible_lamina_id]
                    ## Delta mstruct of the emerged lamina
                    delta_lamina_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['SSLW'], curr_visible_lamina_inputs['mstruct'], curr_visible_lamina_inputs['green_area'])
                    ## Delta Nstruct of the emerged lamina
                    delta_lamina_Nstruct = model.calculate_delta_Nstruct(delta_lamina_mstruct)
                    ## Export of metabolite from hiddenzone towards emerged lamina
                    export_sucrose = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                    export_amino_acids = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])
                    export_fructan = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['fructan'], hiddenzone_inputs['mstruct'])
                    export_proteins = model.calculate_export(delta_lamina_mstruct, hiddenzone_inputs['proteins'], hiddenzone_inputs['mstruct'])

                    # Update of lamina outputs
                    curr_visible_lamina_outputs['mstruct'] += delta_lamina_mstruct
                    curr_visible_lamina_outputs['Nstruct'] += delta_lamina_Nstruct
                    curr_visible_lamina_outputs['sucrose'] += export_sucrose
                    curr_visible_lamina_outputs['amino_acids'] += export_amino_acids
                    curr_visible_lamina_outputs['fructan'] += export_fructan
                    curr_visible_lamina_outputs['proteins'] += export_proteins

                    self.outputs['elements'][visible_lamina_id] = curr_visible_lamina_outputs

                else: #: Mature lamina, growing sheath
                    visible_sheath_id = hiddenzone_id + tuple(['sheath','StemElement']) # The hidden part of the sheath is only updated once, at the end of leaf elongation, by remobilisation from the hiddenzone
                    curr_visible_sheath_inputs = all_elements_inputs[visible_sheath_id]
                    curr_visible_sheath_outputs = all_elements_outputs[visible_sheath_id]
                    ## Delta mstruct of the emerged sheath
                    delta_sheath_mstruct = model.calculate_delta_emerged_tissue_mstruct(hiddenzone_inputs['LSSW'], curr_visible_sheath_inputs['mstruct'], curr_visible_sheath_inputs['length'])
                    ## Delta Nstruct of the emerged sheath
                    delta_sheath_Nstruct = model.calculate_delta_Nstruct(delta_sheath_mstruct)
                    ## Export of metabolite from hiddenzone towards emerged sheath
                    export_sucrose = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['sucrose'], hiddenzone_inputs['mstruct'])
                    export_amino_acids = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['amino_acids'], hiddenzone_inputs['mstruct'])
                    export_fructan = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['fructan'], hiddenzone_inputs['mstruct'])
                    export_proteins = model.calculate_export(delta_sheath_mstruct, hiddenzone_inputs['proteins'], hiddenzone_inputs['mstruct'])

                    # Update of sheath outputs
                    curr_visible_sheath_outputs['mstruct'] += delta_sheath_mstruct
                    curr_visible_sheath_outputs['Nstruct'] += delta_sheath_Nstruct
                    curr_visible_sheath_outputs['sucrose'] += export_sucrose
                    curr_visible_sheath_outputs['amino_acids'] += export_amino_acids
                    curr_visible_sheath_outputs['fructan'] += export_fructan
                    curr_visible_sheath_outputs['proteins'] += export_proteins
                    self.outputs['elements'][visible_sheath_id] = curr_visible_sheath_outputs


            # CN consumption due to mstruct/Nstruct growth of the enclosed leaf and of the internode
            curr_hiddenzone_outputs['AA_consumption_mstruct'] = model.calculate_s_Nstruct_amino_acids( (delta_leaf_enclosed_Nstruct + delta_internode_enclosed_Nstruct), delta_lamina_Nstruct, delta_sheath_Nstruct)  #: Consumption of amino acids due to mstruct growth (µmol N)
            curr_hiddenzone_outputs['sucrose_consumption_mstruct'] = model.calculate_s_mstruct_sucrose( (delta_leaf_enclosed_mstruct + delta_internode_enclosed_mstruct), delta_lamina_mstruct, delta_sheath_mstruct, curr_hiddenzone_outputs['AA_consumption_mstruct']) #: Consumption of sucrose due to mstruct growth (µmol C)
            curr_hiddenzone_outputs['Respi_growth'] = RespirationModel.R_growth(curr_hiddenzone_outputs['sucrose_consumption_mstruct'])                                      #: Respiration growth (µmol C)

            # Update of outputs
            curr_hiddenzone_outputs['leaf_enclosed_mstruct'] += delta_leaf_enclosed_mstruct
            curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] += delta_leaf_enclosed_Nstruct
            curr_hiddenzone_outputs['internode_enclosed_mstruct'] += delta_internode_enclosed_mstruct
            curr_hiddenzone_outputs['internode_enclosed_Nstruct'] += delta_internode_enclosed_Nstruct
            curr_hiddenzone_outputs['mstruct'] = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] + curr_hiddenzone_outputs['internode_enclosed_mstruct']
            curr_hiddenzone_outputs['sucrose'] -= (curr_hiddenzone_outputs['sucrose_consumption_mstruct'] + curr_hiddenzone_outputs['Respi_growth'] + export_sucrose) # TODO: Add checks for negative sucrose value
            curr_hiddenzone_outputs['amino_acids'] -= (curr_hiddenzone_outputs['AA_consumption_mstruct'] + export_amino_acids) # TODO: Add checks for negative AA value
            self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs

            # Remobilisation at the end of leaf elongation
            if not hiddenzone_inputs['leaf_is_growing'] and hiddenzone_inputs['delta_leaf_L'] > 0:
               share = curr_hiddenzone_outputs['leaf_enclosed_mstruct'] / curr_hiddenzone_outputs['mstruct']

               ## Add to hidden part of the sheath
               hidden_sheath_id = hiddenzone_id + tuple(['sheath','HiddenElement'])
               curr_hidden_sheath_outputs = self.outputs['elements'][hidden_sheath_id]
               curr_hidden_sheath_outputs['mstruct'] += curr_hiddenzone_outputs['leaf_enclosed_mstruct']
               curr_hidden_sheath_outputs['Nstruct'] += curr_hiddenzone_outputs['leaf_enclosed_Nstruct']
               curr_hidden_sheath_outputs['sucrose'] += curr_hiddenzone_outputs['sucrose'] * share
               curr_hidden_sheath_outputs['amino_acids'] += curr_hiddenzone_outputs['amino_acids'] * share
               curr_hidden_sheath_outputs['fructan'] += curr_hiddenzone_outputs['fructan'] * share
               curr_hidden_sheath_outputs['proteins'] += curr_hiddenzone_outputs['proteins'] * share
               self.outputs['elements'][hidden_sheath_id] = curr_hidden_sheath_outputs

               ## Remove in hiddenzone
               curr_hiddenzone_outputs = self.outputs['hiddenzone'][hiddenzone_id]
               curr_hiddenzone_outputs['leaf_enclosed_mstruct'] = 0
               curr_hiddenzone_outputs['leaf_enclosed_Nstruct'] = 0
               curr_hiddenzone_outputs['mstruct'] = curr_hiddenzone_outputs['internode_enclosed_mstruct']
               curr_hiddenzone_outputs['Nstruct'] = curr_hiddenzone_outputs['internode_enclosed_Nstruct']
               curr_hiddenzone_outputs['sucrose'] -= curr_hiddenzone_outputs['sucrose'] * share
               curr_hiddenzone_outputs['amino_acids'] -= curr_hiddenzone_outputs['amino_acids'] * share
               curr_hiddenzone_outputs['fructan'] -= curr_hiddenzone_outputs['fructan'] * share
               curr_hiddenzone_outputs['proteins'] -= curr_hiddenzone_outputs['proteins'] * share
               self.outputs['hiddenzone'][hiddenzone_id] = curr_hiddenzone_outputs

            # Remobilisation at the end of internode elongation
            if not hiddenzone_inputs['internode_is_growing'] and hiddenzone_inputs['internode_L'] > 0: # Internodes stop to elongate after leaves. We cannot test delta_internode_L > 0 for the cases of short internodes which are mature before GA production.

                ## Add to hidden part of the internode
                hidden_internode_id = hiddenzone_id + tuple(['internode','HiddenElement'])
                if hidden_internode_id not in self.outputs['elements'].keys():
                    new_internode_outputs = parameters.OrganInit().__dict__
                    self.outputs['elements'][hidden_internode_id] = new_internode_outputs
                curr_hidden_internode_outputs = self.outputs['elements'][hidden_internode_id] # TODO: Initialise internode element outputs when starting its elongation
                curr_hidden_internode_outputs['mstruct'] += curr_hiddenzone_outputs['internode_enclosed_mstruct']
                curr_hidden_internode_outputs['Nstruct'] += curr_hiddenzone_outputs['internode_enclosed_Nstruct']
                curr_hidden_internode_outputs['sucrose'] += curr_hiddenzone_outputs['sucrose']
                curr_hidden_internode_outputs['amino_acids'] += curr_hiddenzone_outputs['amino_acids']
                curr_hidden_internode_outputs['fructan'] += curr_hiddenzone_outputs['fructan']
                curr_hidden_internode_outputs['proteins'] += curr_hiddenzone_outputs['proteins']
                curr_hidden_internode_outputs['is_growing'] = False
                self.outputs['elements'][hidden_internode_id] = curr_hidden_internode_outputs

            #: Delete Hiddenzone after remobilisation so it is not sent to CN Wheat
            if not hiddenzone_inputs['leaf_is_growing'] and not hiddenzone_inputs['internode_is_growing']:
                del self.outputs['hiddenzone'][hiddenzone_id]


        # Roots
        for root_id, root_inputs in all_roots_inputs.iteritems():
            curr_root_outputs = all_roots_outputs[root_id]
            # Growth
            mstruct_C_growth, mstruct_growth, Nstruct_growth, Nstruct_N_growth = model.calculate_roots_mstruct_growth(root_inputs['sucrose'], root_inputs['amino_acids'], root_inputs['mstruct'], self.delta_t)
            # Respiration growth
            curr_root_outputs['Respi_growth'] = RespirationModel.R_growth(mstruct_C_growth)
            # Update of root outputs
            curr_root_outputs['mstruct'] += mstruct_growth
            curr_root_outputs['sucrose'] -= (mstruct_C_growth + curr_root_outputs['Respi_growth'])
            curr_root_outputs['Nstruct'] += Nstruct_growth
            curr_root_outputs['amino_acids'] -= Nstruct_N_growth
            curr_root_outputs['rate_mstruct_growth'] = mstruct_growth / self.delta_t
            self.outputs['roots'][root_id] = curr_root_outputs