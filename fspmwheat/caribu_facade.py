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
        self._c_scene = CaribuScene(pattern=geometrical_model.domain)
        self._shared_elements_inputs_outputs_df = shared_elements_inputs_outputs_df #: the dataframe at elements scale shared between all models
        # Diffuse light sources
        sky = GenSky.GenSky()(1, 'soc', 4, 5) # (Energy, soc/uoc, azimuts, zenits)
        sky = GetLight.GetLight(sky)
        self._c_scene.addSources(sky)
    
    
    def run(self):
        """
        Run the model and update the MTG and the dataframes shared between all models.
        """
        idmap = self._initialize_model()
        output = self._c_scene.runCaribu(direct=True, infinity=True)
        # # Visualisation
        # c_scene.plot(output=output)
        # Aggregation of results
        res_by_id = self._c_scene.output_by_id(output, idmap)
        self._update_shared_MTG(res_by_id)
        self._update_shared_dataframes(res_by_id)
        
        
    def _initialize_model(self):
        """
        Initialize the inputs of the model from the MTG shared between all models and the soils.
        """
        self._c_scene.resetScene()
        geom = self._shared_mtg.property('geometry')
        # Labels
        labels = []
        for vid in geom.keys():
            if self._shared_mtg.class_name(vid) == 'HiddenElement':
                continue
            elif self._shared_mtg.class_name(vid) in ('LeafElement1', 'LeafElement'):
                plant_id = self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, scale=1))
                label = encode_label(opt_id=1, opak=1, plant_id=plant_id)[0]
            elif self._shared_mtg.class_name(vid) == 'StemElement':
                plant_id = self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, scale=1))
                label = encode_label(opt_id=1, opak=0, plant_id=plant_id)[0]
            else:
                label = encode_label()[0]
                print 'Warning: unknown element type {}, vid={}'.format(self._shared_mtg.class_name(vid), vid)
            labels.append(label)

        # Add scene to CaribuScene
        idmap = self._c_scene.add_Shapes(geom, canlabels=labels)
        return idmap
    
    
    def _update_shared_MTG(self, caribu_results_by_id):
        """
        Update the MTG shared between all models from the population of Caribu.
        """
        # add the missing property
        if 'Eabsm2' not in self._shared_mtg.properties():
            self._shared_mtg.add_property('Eabsm2')
        
        # update the MTG
        self._shared_mtg.property('Eabsm2').update(caribu_results_by_id['Eabsm2'])
        
    
    def _update_shared_dataframes(self, caribu_results_by_id):
        """
        Update the dataframes shared between all models from the inputs dataframes or the outputs dataframes of the model.
        """
        ids = []
        for vid in caribu_results_by_id['label'].keys():
            ind = (int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 1))), 
                   self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 2)), 
                   int(self._shared_mtg.index(self._shared_mtg.complex_at_scale(vid, 3))), 
                   self._shared_mtg.label(self._shared_mtg.complex_at_scale(vid, 4)), 
                   self._shared_mtg.label(vid))
            ids.append(ind)
    
        ids_df = pd.DataFrame(sorted(ids), columns=SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES)
        Eabsm2_list = np.array(sorted(caribu_results_by_id['Eabsm2'].items()))[:,1] # Used in order to sort values by keys.
        data_df = pd.DataFrame({'Eabsm2': Eabsm2_list})
        df = pd.concat([ids_df, data_df], axis=1)
        df.reset_index(drop=True, inplace=True)
        tools.combine_dataframes_inplace(df, SHARED_ELEMENTS_INPUTS_OUTPUTS_INDEXES, self._shared_elements_inputs_outputs_df)
        
