# -*- coding: latin-1 -*-
"""
    fspmwheat
    ~~~~~~~

    The model FSPM-Wheat.

    A Functional Structural Plant Model of Wheat.

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

__version__  = '0.0.1'

# Add a do-nothing handler to prevent an error message being output to sys.stderr in the absence of logging configuration
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())