##########################################################################
# NSAp - Copyright (C) CEA, 2019
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


"""
Module for generating post-preproc plots (registration, segmentation, etc.)
"""


# Imports
import numpy as np
import matplotlib.pyplot as plt


def plot_fsl_motion_parameters(parameter_file, outfname):
    """ Plot motion parameters obtained with FSL software

    Parameters
    ----------
    parameter_file: string
        path of file containing the motion parameters.
    outfname: string
        output filename for storing the plotted figure.

    """
    # Load parameters
    motion = np.loadtxt(parameter_file)
    motion[:, :3] *= (180. / np.pi)

    # do plotting
    plt.figure()
    fig, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(motion[:, 3:])
    ax1.set_xlabel("time(scans)")
    ax1.set_ylabel("Estimated motion (mm)")
    ax1.grid(True)
    ax2.plot(motion[:, :3])
    ax2.set_xlabel("time(scans)")
    ax2.set_ylabel("Estimated motion (degrees)")
    ax2.grid(True)
    plt.legend(("TransX", "TransY", "TransZ", "RotX", "RotY", "RotZ"),
               loc="upper left", ncol=2)
    plt.savefig(outfname)
    plt.close()

