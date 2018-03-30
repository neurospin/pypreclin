##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Package to process fMRI datasets.
"""

from .info import __version__
from .info import DEFAULT_SPM_STANDALONE_PATH
from .info import DEFAULT_FSL_PATH
from .configure import info


# Display a welcome message
print(info())
