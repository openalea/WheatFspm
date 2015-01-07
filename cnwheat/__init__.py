# -*- coding: latin-1 -*-
"""
    cnwheat
    ~~~~~~~

    The model CN-Wheat.
    
    CN-Wheat computes CN exchanges in wheat architecture defined by laminae, 
    sheaths, internodes, a peduncle, a chaff, a phloem, roots and grains.

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

__version__  = '0.0.1'

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())