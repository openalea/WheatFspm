# -*- coding: latin-1 -*-

from __future__ import division  # use "//" to do integer division

import model

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


class SimulationError(Exception): pass


class SimulationInputsError(SimulationError): pass


class Simulation(object):
    """The Simulation class permits to initialize and run a simulation.
    """

    def __init__(self):

        #: The inputs of Farquhar-Wheat.
        #:
        #: `inputs` is a dictionary of dictionaries:
        #:     {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_input_name: element_input_value, ...}, ...}
        #: See :meth:`Model.run <farquharwheat.model.Model.run>`
        #: for more information about the inputs.
        self.inputs = {}

        #: The outputs of Farquhar-Wheat.
        #:
        #: `outputs` is a dictionary of dictionaries:
        #:     {(plant_index, axis_label, metamer_index, organ_label, element_label): {element_output_name: element_output_value, ...}, ...}
        #: See :meth:`Model.run <farquharwheat.model.Model.run>`
        #: for more information about the outputs.
        self.outputs = {}

    def initialize(self, inputs):
        """
        Initialize :attr:`inputs` from `inputs`.

        :Parameters:

            - `inputs` (:class:`dict`) - Dictionary of two dictionaries :
                    - `elements` : The inputs by element.
                    - `SAMs` : The inputs by axis.
              `inputs` must be a dictionary with the same structure as :attr:`inputs`.

            See :meth:`Model.run <farquharwheat.model.Model.run>`
               for more information about the inputs.
        """
        self.inputs.clear()
        self.inputs.update(inputs)

    def run(self, Ta, ambient_CO2, RH, Ur):
        """
        Compute Farquhar variables for each element in :attr:`inputs` and put
        the results in :attr:`outputs`.

        :Parameters:

            - `Ta` (:class:`float`) - air temperature at t (degree Celsius)

            - `ambient_CO2` (:class:`float`) - air CO2 at t (µmol mol-1)

            - `RH` (:class:`float`) - relative humidity at t (decimal fraction)

            - `Ur` (:class:`float`) - wind speed at the top of the canopy at t (m s-1)

        """
        self.outputs.update({inputs_type: {} for inputs_type in self.inputs['elements'].keys()})

        for (element_id, element_inputs) in self.inputs['elements'].items():

            SAM_id = element_id[:2]
            organ_label = element_id[3]
            element_label = element_id[4]

            axe_label = SAM_id[1]
            if axe_label != 'MS': # Calculation only for the main stem
                continue

            if element_inputs['height'] is None or element_inputs['green_area'] == 0.0:  # In case it is an HiddenElement, we need temperature calculation. Cases of Visible Element without geomtry proprety (because too small) don't have photosynthesis calculation neither.
                #element_label == 'HiddenElement' or
                Ag, An, Rd, Tr, gs = 0.0, 0.0, 0.0, 0.0, 0.0
                Ts = self.inputs['SAMs'][SAM_id]['SAM_temperature']
                Tr = 0.1 # Default transpiration value for small organs under ADEL's resolution (green_area == 0)
            else:
                PARa = element_inputs['PARa']     #: Amount of absorbed PAR per unit area (µmol m-2 s-1)

                surfacic_nitrogen = model.Model.calculate_surfacic_nitrogen(element_inputs['nitrates'],
                                                                        element_inputs['amino_acids'],
                                                                        element_inputs['proteins'],
                                                                        element_inputs['Nstruct'],
                                                                        element_inputs['green_area'])

                height_canopy = self.inputs['SAMs'][SAM_id]['height_canopy']
                Ag, An, Rd, Tr, Ts, gs = model.Model.run(surfacic_nitrogen,
                                                                  element_inputs['width'],
                                                                  element_inputs['height'],
                                                                  PARa, Ta, ambient_CO2, RH, Ur, organ_label, height_canopy)

            element_outputs = {'Ag': Ag, 'An': An, 'Rd': Rd, 'Tr': Tr, 'Ts': Ts, 'gs': gs, 'width': element_inputs['width'], 'height': element_inputs['height']}

            self.outputs[element_id] = element_outputs
