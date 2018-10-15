# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

import model

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


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    MIN_GREEN_AREA = 0.5E-4  #: Minimal green area of an element (m2). Below this area, set green_area to 0.0.

    def __init__(self, delta_t=1):

        #: The inputs of Senesc-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {'roots': {(plant_index, axis_label): {roots_input_name: roots_input_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_input_name: element_input_value, ...}, ...}}
        self.inputs = {}

        #: The outputs of Senesc-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {'roots': {(plant_index, axis_label): {roots_output_name: roots_output_value, ...}, ...},
        #:      'elements': {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_output_name: element_output_value, ...}, ...}}
        self.outputs = {}

        #: the delta t of the simulation (in seconds)
        self.delta_t = delta_t

    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :Parameters:

            - `inputs` (:class:`dict`) - The inputs by roots and element.
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self, forced_max_protein_elements=None):
        """
        Compute Senesc-Wheat outputs from :attr:`inputs`, and update :attr:`outputs`.

        :Parameters:

            - `forced_max_protein_elements` (:class:`set`) - The elements ids with fixed max proteins.

        .. todo:: remove forced_max_protein_elements

        """
        self.outputs.update({inputs_type: {} for inputs_type in self.inputs.keys()})

        # SAM
        all_SAM_inputs = self.inputs['SAM']

        # Roots
        all_roots_inputs = self.inputs['roots']
        all_roots_outputs = self.outputs['roots']
        for roots_inputs_id, roots_inputs_dict in all_roots_inputs.items():

            # Temperature-compensated time (delta_teq)
            delta_teq = all_SAM_inputs[roots_inputs_id]['delta_teq_roots']

            # loss of mstruct and Nstruct
            rate_mstruct_death, rate_Nstruct_death = model.SenescenceModel.calculate_roots_senescence(roots_inputs_dict['mstruct'], roots_inputs_dict['Nstruct'])
            relative_delta_mstruct = model.SenescenceModel.calculate_relative_delta_mstruct_roots(rate_mstruct_death, roots_inputs_dict['mstruct'], delta_teq)
            # loss of cytokinins (losses of nitrates, amino acids and sucrose are neglected)
            loss_cytokinins = model.SenescenceModel.calculate_remobilisation(roots_inputs_dict['cytokinins'], relative_delta_mstruct)
            # Update of root outputs
            all_roots_outputs[roots_inputs_id] = {'mstruct': roots_inputs_dict['mstruct'] - (rate_mstruct_death * delta_teq),  #: TODO: a faire dans une fonction a part et apres growth-wheat
                                                  'rate_mstruct_death': rate_mstruct_death,
                                                  'Nstruct': roots_inputs_dict['Nstruct'] - (rate_Nstruct_death * delta_teq),  #: TODO: a faire dans une fonction a part et apres growth-wheat
                                                  'cytokinins': roots_inputs_dict['cytokinins'] - loss_cytokinins}


        # Elements
        all_elements_inputs = self.inputs['elements']
        all_elements_outputs = self.outputs['elements']
        for element_inputs_id, element_inputs_dict in all_elements_inputs.items():

            axe_label = element_inputs_id[1]
            if axe_label != 'MS': # Calculation only for the main stem
                continue

            # Temperature-compensated time (delta_teq)
            axe_id = element_inputs_id[:2]
            delta_teq = all_SAM_inputs[axe_id]['delta_teq']

            # Senescence
            if element_inputs_dict['green_area'] < Simulation.MIN_GREEN_AREA and not element_inputs_dict['is_growing']:
                element_outputs_dict = element_inputs_dict.copy()
                element_outputs_dict['green_area'] = 0.0
            else:
                update_max_protein = forced_max_protein_elements is None or not element_inputs_id in forced_max_protein_elements
                # new_green_area, relative_delta_green_area, max_proteins = model.SenescenceModel.calculate_relative_delta_green_area(element_inputs_id[3], element_inputs_dict['green_area'],
                #                                                                                                                     element_inputs_dict['proteins'] / element_inputs_dict['mstruct'],
                #                                                                                                                     element_inputs_dict['max_proteins'], self.delta_t, update_max_protein)
                # Temporaire
                new_senesced_length, relative_delta_senesced_length, max_proteins = model.SenescenceModel.calculate_relative_delta_senesced_length(element_inputs_id[3],
                                                                                                                                                element_inputs_dict['senesced_length'],
                                                                                                                                                element_inputs_dict['length'],
                                                                                                                                                element_inputs_dict['proteins'] / element_inputs_dict['mstruct'],
                                                                                                                                                element_inputs_dict['max_proteins'], delta_teq,
                                                                                                                                                update_max_protein)
                # Temporaire :
                relative_delta_green_area = relative_delta_senesced_length

                new_mstruct, new_Nstruct = model.SenescenceModel.calculate_delta_mstruct_shoot(relative_delta_green_area, element_inputs_dict['mstruct'], element_inputs_dict['Nstruct'])

                # Remobilisation
                remob_starch = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['starch'], relative_delta_green_area)
                remob_fructan = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['fructan'], relative_delta_green_area)
                remob_proteins = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['proteins'], relative_delta_green_area)
                loss_cytokinins = model.SenescenceModel.calculate_remobilisation(element_inputs_dict['cytokinins'], relative_delta_green_area)

                element_outputs_dict = {'green_area': element_inputs_dict['green_area'],#new_green_area,
                                        'senesced_length' : new_senesced_length,
                                        'mstruct': new_mstruct,
                                        'Nstruct': new_Nstruct,
                                        'starch': element_inputs_dict['starch'] - remob_starch,
                                        'sucrose': element_inputs_dict['sucrose'] + remob_starch + remob_fructan,
                                        'fructan': element_inputs_dict['fructan'] - remob_fructan,
                                        'proteins': element_inputs_dict['proteins'] - remob_proteins,
                                        'amino_acids': element_inputs_dict['amino_acids'] + remob_proteins,
                                        'cytokinins': element_inputs_dict['cytokinins'] - loss_cytokinins,
                                        'max_proteins': max_proteins}

            all_elements_outputs[element_inputs_id] = element_outputs_dict
