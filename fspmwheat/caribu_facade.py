# -*- coding: latin-1 -*-

"""
    fspmwheat.caribu_facade
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`fspmwheat.caribu_facade` is a facade of the model Caribu.

    This module permits to initialize and run the model Caribu from a :class:`MTG <openalea.mtg.mtg.MTG>`
    in a convenient and transparent way, wrapping all the internal complexity of the model, and dealing
    with all the tedious initialization and conversion processes.

    :copyright: Copyright 2014-2016 INRA-ECOSYS, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2016.
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

from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.sky_tools import GenSky, GetLight
from alinea.caribu.label import encode_label

import tools

#: the columns which define the topology in the elements scale dataframe shared between all models
SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES = ['plant', 'axis', 'metamer', 'organ', 'element']

#: the outputs of Caribu
CARIBU_OUTPUTS = ['Eabsm2']


class CaribuFacade(object):
    """
    The CaribuFacade class permits to initialize, run the model Caribu
    from a :class:`MTG <openalea.mtg.mtg.MTG>`, and update the MTG and the dataframes
    shared between all models.

    Use :meth:`run` to run the model.

    :Parameters:

        - `shared_mtg` (:class:`openalea.mtg.mtg.MTG`) - The MTG shared between all models.
        - `shared_elements_inputs_outputs_df` (:class:`pandas.DataFrame`) - the dataframe of inputs and outputs at elements scale shared between all models.
        - `geometrical_model` (:func:`geometrical_model`) - The model which deals with geometry. This model must have an attribute "domain".

    """

    def __init__(self,
                 shared_mtg,
                 shared_elements_inputs_outputs_df,
                 geometrical_model):

        self._shared_mtg = shared_mtg #: the MTG shared between all models
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df #: the dataframe at elements scale shared between all models
        self._geometrical_model = geometrical_model #: the model which deals with geometry

    def run(self):
        """
        Run the model and update the MTG and the dataframes shared between all models.
        """
        c_scene = self._initialize_model()
        _, aggregated = c_scene.run(direct=True, infinite=True)

##        # Visualisation
##        c_scene.plot()

        # Eabs is the absorbed energy per m2
        self._update_shared_MTG(aggregated['par']['Eabs'])
        self._update_shared_dataframes(aggregated['par']['Eabs'])

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models and the soils.
        """

        # Diffuse light sources
        sky = GenSky.GenSky()(1, 'soc', 4, 5) # (Energy, soc/uoc, azimuts, zenits)
        sky = GetLight.GetLight(sky)

        # Optical
        opt = {'par':{}}
        geom = self._shared_mtg.property('geometry')
        for vid in geom.keys():
            if self._shared_mtg.class_name(vid) == 'HiddenElement':#: TODO: check if necessary?
                continue
            elif self._shared_mtg.class_name(vid) in ('LeafElement1', 'LeafElement'):
                opt['par'][vid] = (0.10, 0.05)
            elif self._shared_mtg.class_name(vid) == 'StemElement':
                opt['par'][vid] = (0.10,)
            else:
                print 'Warning: unknown element type {}, vid={}'.format(self._shared_mtg.class_name(vid), vid)

        c_scene = CaribuScene(scene=self._shared_mtg, light=sky, pattern=self._geometrical_model.domain, opt=opt)

        return c_scene

    def _update_shared_MTG(self, aggregated_Eabsm2):
        """
        Update the MTG shared between all models from the population of Caribu.
        """
        # add the missing property
        if 'Eabsm2' not in self._shared_mtg.properties():
            self._shared_mtg.add_property('Eabsm2')

        # update the MTG
        self._shared_mtg.property('Eabsm2').update(aggregated_Eabsm2)


    def _update_shared_dataframes(self, aggregated_Eabsm2):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        """
        ids = []
        Eabsm2 = []
        for vid in sorted(aggregated_Eabsm2.keys()):
            if int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 1))) ==1 and self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 2)) == 'MS': #TODO: temporary patch for test
                ind = int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 1))), self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 2)), int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 3))), self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 4)), self._shared_mtg.label(vid)
                ids.append(ind)
                Eabsm2.append(aggregated_Eabsm2[vid])

        ids_df = pd.DataFrame(ids, columns=SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES)
        data_df = pd.DataFrame({'Eabsm2': Eabsm2})
        df = pd.concat([ids_df, data_df], axis=1)
        df.sort_values(['plant', 'axis', 'metamer', 'organ', 'element'], inplace=True)
        tools.combine_dataframes_inplace(df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)

