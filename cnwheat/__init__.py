# -*- coding: latin-1 -*-
"""
    cnwheat
    ~~~~~~~

    The model CN-Wheat.
    
    CN-Wheat computes the CN exchanges in a wheat architecture. See:
    
        * :mod:`cnwheat.simulation` for the front-end of the model,
        * :mod:`cnwheat.model` for the equations of the model,
        * :mod:`cnwheat.parameters` for the parameters of the model.
        
    :copyright: Copyright 2014-2017 INRA-ECOSYS, see AUTHORS.
    :license: CeCILL-C, see LICENSE for details.
    
    **Acknowledgments**: The research leading these results has received funding through the 
    Investment for the Future programme managed by the Research National Agency 
    (BreedWheat project ANR-10-BTBR-03).
    
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

__version__  = '0.0.1'

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())