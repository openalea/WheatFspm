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

    def run(self, PARi):
        """
        Run the model and update the MTG and the dataframes shared between all models.

        :Parameters:
        - `PARi` (:class:`float`) - Incident PAR above the canopy (µmol m-2 s-1)

        """
        c_scene = self._initialize_model()
        _, aggregated = c_scene.run(direct=True, infinite=True)

        Eabs = aggregated['par']['Eabs']
        Eabs.update((vid, Eabs * PARi) for vid, Eabs in Eabs.iteritems())

##        # Visualisation
##        c_scene.plot()

        # Eabs is the relative surfacic absorbed energy per organ
        self._update_shared_MTG(aggregated['par']['Eabs'])
        self._update_shared_dataframes(aggregated['par']['Eabs'])

    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared
        """

        # Diffuse light sources
        sky_string = GetLight.GetLight(GenSky.GenSky()(1, 'soc', 4, 5)) # (Energy, soc/uoc, azimuts, zenits)

        sky = []
        for string in sky_string.split('\n'):
            if len(string)!=0:
                string_split = string.split(' ')
                t = tuple((float(string_split[0]), tuple((float(string_split[1]), float(string_split[2]), float(string_split[3])))))
                sky.append(t)

        # Optical
        opt = {'par':{}}
        geom = self._shared_mtg.property('geometry')

        for vid in geom.keys():
            if self._shared_mtg.class_name(vid) in ('LeafElement1', 'LeafElement'):
                opt['par'][vid] = (0.10, 0.05) #: (reflectance, transmittance) of the adaxial side of the leaves
            elif self._shared_mtg.class_name(vid) == 'StemElement':
                opt['par'][vid] = (0.10,) #: (reflectance,) of the stems
            else:
                print 'Warning: unknown element type {}, vid={}'.format(self._shared_mtg.class_name(vid), vid)

        c_scene = CaribuScene(scene=self._shared_mtg, light=sky, pattern=self._geometrical_model.domain, opt=opt)

        return c_scene

    def _update_shared_MTG(self, aggregated_PARa):
        """
        Update the MTG shared between all models from the population of Caribu.
        """
        # add the missing property
        if 'PARa' not in self._shared_mtg.properties():
            self._shared_mtg.add_property('PARa')

        # update the MTG
        self._shared_mtg.property('PARa').update(aggregated_PARa)


    def _update_shared_dataframes(self, aggregated_PARa):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        """
        ids = []
        PARa = []
        for vid in sorted(aggregated_PARa.keys()):
            ind = int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 1))), self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 2)), int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 3))), self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 4)), self._shared_mtg.label(vid)
            ids.append(ind)
            PARa.append(aggregated_PARa[vid])

        ids_df = pd.DataFrame(ids, columns=SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES)
        data_df = pd.DataFrame({'PARa': PARa})
        df = pd.concat([ids_df, data_df], axis=1)
        df.sort_values(['plant', 'axis', 'metamer', 'organ', 'element'], inplace=True)
        tools.combine_dataframes_inplace(df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)

    def run_from_df(self, Eabs_df, PARi, multiple_sources=False, ratio_diffus_PAR=None):
        """
        Update the MTG and the dataframes shared between all models from an input dataframe having Eabms2 values.
        """


        if 'species' in self._shared_mtg.properties():
            Eabs_df_grouped = Eabs_df.groupby(['species', 'metamer', 'organ'])
        else:
            Eabs_df_grouped = Eabs_df.groupby(['metamer', 'organ'])

        #: the name of the organs modeled by FarquharWheat
        CARIBU_ORGANS_NAMES = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])
        #: the name of the elements modeled by FarquharWheat
        CARIBU_ELEMENTS_NAMES = set(['StemElement', 'LeafElement1'])

        PARa_element_data_dict = {}
        # traverse the MTG recursively from top ...
        for mtg_plant_vid in self._shared_mtg.components_iter(self._shared_mtg.root):
            for mtg_axis_vid in self._shared_mtg.components_iter(mtg_plant_vid):
                for mtg_metamer_vid in self._shared_mtg.components_iter(mtg_axis_vid):
                    mtg_metamer_index = int(self._shared_mtg.index(mtg_metamer_vid))
                    for mtg_organ_vid in self._shared_mtg.components_iter(mtg_metamer_vid):
                        mtg_organ_label = self._shared_mtg.label(mtg_organ_vid)
                        for mtg_element_vid in self._shared_mtg.components_iter(mtg_organ_vid):
                            mtg_element_label = self._shared_mtg.label(mtg_element_vid)
                            if mtg_element_label not in CARIBU_ELEMENTS_NAMES: continue
                            if 'species' in self._shared_mtg.properties():
                                element_id = (self._shared_mtg.property('species')[mtg_plant_vid], mtg_metamer_index, mtg_organ_label)
                            else:
                                element_id = (mtg_metamer_index, mtg_organ_label)
                            if element_id in Eabs_df_grouped.groups.keys():
                                if PARi == 0:
                                    PARa_element_data_dict[mtg_element_vid] = 0
                                elif multiple_sources:
                                    PARa_diffuse = Eabs_df_grouped.get_group(element_id)['Eabs_diffuse'].iloc[0] * PARi * ratio_diffus_PAR
                                    PARa_direct = Eabs_df_grouped.get_group(element_id)['Eabs_direct'].iloc[0] * PARi * (1 - ratio_diffus_PAR)
                                    PARa_element_data_dict[mtg_element_vid] = PARa_diffuse + PARa_direct
                                else:
                                    PARa_element_data_dict[mtg_element_vid] = Eabs_df_grouped.get_group(element_id)['Eabs'].iloc[0] * PARi

        # update MTG and datagrame
        self._update_shared_MTG(PARa_element_data_dict)
        self._update_shared_dataframes(PARa_element_data_dict)