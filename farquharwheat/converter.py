# -*- coding: latin-1 -*-

from __future__ import division # use "//" to do integer division

"""
    farquharwheat.converter
    ~~~~~~~~~~~~~~~~~~~~~~~

    The module :mod:`farquharwheat.converter` permits to:
        
        * convert a :class:`MTG <openalea.mtg.mtg.MTG>` to Farquhar-Wheat inputs dictionary (see :func:`from_MTG`),
        * update a :class:`MTG <openalea.mtg.mtg.MTG>` from Farquhar-Wheat outputs dictionary (see :func:`update_MTG`).

    :copyright: Copyright 2014 INRA-EGC, see AUTHORS.
    :license: TODO, see LICENSE for details.

    .. seealso:: Barillot et al. 2014.
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

from alinea.astk import plantgl_utils

import warnings

class ConverterWarning(UserWarning): pass

class PropertyNotFoundWarning(ConverterWarning):
    '''Property not found in a vertex of a MTG.'''
    def __init__(self, property_, vertex_id):
        self.message = 'Property {0} not found in vertex {1}. Ignore vertex {1} and all its components.'.format(property_, vertex_id)
    def __str__(self):
        return repr(self.message)

class InvalidMTGWarning(ConverterWarning):
    '''The input MTG does not contain the required properties.'''
    pass

warnings.simplefilter('always', ConverterWarning)

#: the reduction factor to apply on the STAR of the visible sheath 
#: element of the precedent metamer to obtain the STAR of the hidden 
#: elements of the current metamer. 
SHEATH_STAR_TRANSMISSION = 0.1 

#: the name of the organs modeled by FarquharWheat
ORGANS_NAMES_SET = set(['internode', 'blade', 'sheath', 'peduncle', 'ear'])

#: the name of the inputs FarquharWheat needs to run
FARQUHARWHEAT_MANDATORY_INPUTS_NAMES = set(['width', 'diameter', 'geometry', 'exposed_area', 'area', 'surfacic_nitrogen'])

#: the name of the outputs of FarquharWheat ; :func:`update_MTG` adds/updates only these outputs. 
FARQUHARWHEAT_OUTPUTS_NAMES = set(['Ag', 'An', 'Rd', 'Tr', 'Torg', 'gs'])


def from_MTG(g):
    """
    Extract the inputs of Farquhar-Wheat from a MTG. 
    
    :Parameters:
        
            - g (:class:`openalea.mtg.mtg.MTG`) - A MTG which contains the inputs
              needed by Farquhar-Wheat to be run on each element.
              
    :Returns:
        The inputs of FarquharWheat by element.
        
        The inputs is a dictionary of dictionaries: {element1_id: element1_inputs, element2_id: element2_inputs, ..., elementN_id: elementN_inputs}, where:
         
            * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_type, element_type),
            * and elementi_inputs is a dictionary: {'elementi_input1_name': elementi_input1_value, 'elementi_input2_name': elementi_input2_value, ..., 'elementi_inputN_name': elementi_inputN_value}.
         
        See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
        for more information about the inputs. 
        
    :Returns Type:
        :class:`dict` of :class:`dict`
        
    .. todo:: we currently use alinea.astk.plantgl_utils.get_height to compute the height of each element. 
        This is time consuming and adds a strong dependency. 
        TODO: think about another way to compute the height of the elements.
              
    """
    elements_inputs = {}
    
    # check needed properties
    if not FARQUHARWHEAT_MANDATORY_INPUTS_NAMES.issubset(g.properties()):
        warnings.warn('The input MTG does not contain the required properties ({}): ignore all the MTG.'.format(FARQUHARWHEAT_MANDATORY_INPUTS_NAMES))
        return elements_inputs
    
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            precedent_metamer_vid = None
            precedent_sheath_has_visible = False
            precedent_sheath_STAR = None
            precedent_sheath_height = None
            axis_id = g.label(axis_vid)
            for metamer_vid in g.components_iter(axis_vid):
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in ORGANS_NAMES_SET:
                        warnings.warn('FarquharWheat does not model organ {}. Ignore vertex {} and all its components. \
The organs modeled by FarquharWheat are: {}.'.format(organ_label, organ_vid, ORGANS_NAMES_SET))
                        continue
                    organ_type = organ_label # get organ type
                    vertex_properties = g.get_vertex_property(organ_vid)
                    if organ_type == 'blade':
                        if 'width' not in vertex_properties:
                            warnings.warn(PropertyNotFoundWarning('width', organ_vid))
                            continue
                        organ_width = vertex_properties['width'] # get organ width
                    else:
                        if 'diameter' not in vertex_properties:
                            warnings.warn(PropertyNotFoundWarning('diameter', organ_vid))
                            continue
                        organ_width = vertex_properties['diameter'] # get organ width
                    number_of_visible_elements = 0
                    for element_vid in g.components_iter(organ_vid):
                        vertex_properties = g.get_vertex_property(element_vid)
                        element_label = g.label(element_vid)
                        if element_label.startswith('Hidden'):
                            element_type = 'hidden'
                            if not precedent_sheath_has_visible: 
                                warnings.warn('The sheath of the precedent metamer (vid={}) does not have \
any visible element: no STAR nor height to use for the hidden elements of the organs of the current metamer (vid={}). \
Ignore current element (vid={}) and all its components.'.format(precedent_metamer_vid, metamer_vid, element_vid))
                                continue
                            organ_height = precedent_sheath_height
                            STAR = precedent_sheath_STAR * SHEATH_STAR_TRANSMISSION
                        else:
                            element_type = 'visible'
                            number_of_visible_elements += 1
                            if 'geometry' not in vertex_properties:
                                warnings.warn(PropertyNotFoundWarning('geometry', element_vid))
                                continue
                            # compute organ height
                            triangles = vertex_properties['geometry']
                            triangles_heights = plantgl_utils.get_height({element_vid:triangles}).values()
                            organ_height = np.mean(triangles_heights)
                            if 'exposed_area' not in vertex_properties:
                                warnings.warn(PropertyNotFoundWarning('exposed_area', element_vid))
                                continue
                            exposed_area = vertex_properties['exposed_area']
                            if 'area' not in vertex_properties:
                                warnings.warn(PropertyNotFoundWarning('area', element_vid))
                                continue
                            # compute STAR
                            area = vertex_properties['area']
                            STAR = exposed_area / area
                            if organ_type == 'sheath': 
                                # keep the height and the STAR for the hidden elements of the next metamer
                                precedent_sheath_has_visible = True
                                precedent_sheath_height = organ_height
                                precedent_sheath_STAR = STAR
                        if 'surfacic_nitrogen' not in vertex_properties:
                            warnings.warn(PropertyNotFoundWarning('surfacic_nitrogen', element_vid))
                            continue
                        Na = vertex_properties['surfacic_nitrogen'] # get Na
                        elements_inputs[(plant_index, axis_id, metamer_index, organ_type, element_type)] = {'organ': organ_type,
                                                                                                            'STAR': STAR, 
                                                                                                            'Na': Na, 
                                                                                                            'organ_width': organ_width, 
                                                                                                            'organ_height': organ_height}
                    if organ_type == 'sheath':
                        if number_of_visible_elements == 0:
                            precedent_sheath_has_visible = False
                            precedent_sheath_height = None
                            precedent_sheath_STAR = None
                        elif number_of_visible_elements > 1:
                            warnings.warn('The sheath of the current metamer (vid={}) has more than one visible elements. \
The STAR and height of one of these elements will be used for the hidden elements of the organs of the next metamer. Which element \
is used is undetermined.'.format(metamer_vid))
                
                precedent_metamer_vid = metamer_vid
                            
    return elements_inputs


def update_MTG(outputs, g):
    """
    Update a MTG from Farquhar-Wheat outputs. 
    The list of Farquhar-Wheat outputs is :mod:`FARQUHARWHEAT_OUTPUTS_NAMES`.
    
    :Parameters:
        
            - outputs (:class:`dict` of :class:`dict`) - The outputs by element. 
        
            `outputs` is a dictionary of dictionaries: {element1_id: element1_outputs, element2_id: element2_outputs, ..., elementN_id: elementN_outputs}, where:
             
                * elementi_id is a tuple: (plant_index, axis_id, metamer_index, organ_type, element_type),
                * and elementi_outputs is a dictionary: {'elementi_output1_name': elementi_output1_value, 'elementi_output2_name': elementi_output2_value, ..., 'elementi_outputN_name': elementi_outputN_value}.
        
            See :meth:`Model.calculate_An <farquharwheat.model.Model.calculate_An>` 
            for more information about the outputs.
        
            - `g` (:class:`openalea.mtg.mtg.MTG`) - The MTG to update from the `outputs` of FarquharWheat. 

    """
    # add the properties if needed
    properties = g.properties()
    for farquharwheat_output_name in FARQUHARWHEAT_OUTPUTS_NAMES:
        if farquharwheat_output_name not in properties:
            g.add_property(farquharwheat_output_name)
    # traverse the MTG recursively from top
    for plant_vid in g.components_iter(g.root):
        plant_index = int(g.index(plant_vid))
        for axis_vid in g.components_iter(plant_vid):
            axis_id = g.label(axis_vid)
            for metamer_vid in g.components(axis_vid): 
                metamer_index = int(g.index(metamer_vid))
                for organ_vid in g.components_iter(metamer_vid):
                    organ_label = g.label(organ_vid)
                    if organ_label not in ORGANS_NAMES_SET:
                        continue
                    organ_type = organ_label
                    for element_vid in g.components_iter(organ_vid):
                        element_label = g.label(element_vid)
                        if element_label.startswith('Hidden'):
                            element_type = 'hidden'
                        else:
                            element_type = 'visible'
                        element_id = (plant_index, axis_id, metamer_index, organ_type, element_type)
                        if element_id in outputs:
                            element_outputs = outputs[element_id]
                            for output_name, output_value in element_outputs.iteritems():
                                # set the property
                                g.property(output_name)[element_vid] = output_value

