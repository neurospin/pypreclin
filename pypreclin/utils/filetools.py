##########################################################################
# NSAP - Copyright (C) CEA, 2019
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
File utilities.
"""

# Sytem import
import os

# Third party import
from pyconnectome import DEFAULT_FSL_PATH
from pyconnectome.wrapper import FSLWrapper


def average_timeserie(input_file, output_file, fslconfig=DEFAULT_FSL_PATH):
    """ Average a timeserie.

    Parameters
    ----------
    input_file: str
        the image to average.
    output_file: str
        the averaged image.
    fslconfig: str, default DEFAULT_FSL_PATH
        the FSL configuration batch.
    """
    # Check the input parameter
    for path in (input_file, ):
        if not os.path.isfile(path):
            raise ValueError("'{0}' is not a valid input file.".format(path))

    # Define the FSL command
    cmd = ["fslmaths", input_file, "-Tmean", output_file]

    # Call fslmaths
    fslprocess = FSLWrapper(shfile=fslconfig)
    fslprocess(cmd=cmd)

