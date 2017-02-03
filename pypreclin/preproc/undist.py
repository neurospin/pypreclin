##########################################################################
# NSAp - Copyright (C) CEA, 2013 - 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# Sytem import
import os
import shutil
import subprocess

# Pyconnectome import
import pypreclin
from pyconnectome.utils.segtools import bet2

# Package import
from pypreclin.utils.export import ungzip_file


def mc_undist(funcfile, outdir, spmpath, fslpath):
    """ Slice by slice motion correction.
    """
    # Extract the brain with FSL bet
    outputs = bet2(funcfile, outdir, outline=False, mask=True,
                   skull=False, nooutput=True, f=0.15, g=0, radius=None,
                   smooth=None, c=None, threshold=False, mesh=False,
                   shfile=fslpath)
    gzmaskfile = outputs[1].replace("_mask", "_msk")
    os.rename(outputs[1], gzmaskfile)
    maskfile = ungzip_file(gzmaskfile, prefix="")
    os.remove(gzmaskfile)

    # Generate the matlab job
    shutil.copy(funcfile, outdir)
    jobfile = os.path.join(
        os.path.dirname(pypreclin.__file__), "preproc", "resources", "leuven",
                        "mc_undist_wouterdante_tasserie_job.m")
    info = {
        "SPMPATH": spmpath,
        "PREPROCPATH": os.path.join(os.path.dirname(pypreclin.__file__),
                                    "preproc"),
        "CWD": outdir + os.path.sep,
        "BASENAME": os.path.basename(funcfile)
    }
    with open(jobfile, "rt") as open_file:
        job = open_file.read()
    jobfile = os.path.join(outdir, "mc_undist_job.m")
    with open(jobfile, "wt") as open_file:
        open_file.write(job.format(**info))

    # Call the undist script
    cmd = ["matlab", "-nodisplay", "-nosplash", "-nodesktop", "-r",
           "run {0}".format(jobfile)]
    subprocess.check_call(cmd)
    undist_file = os.path.join(outdir, "u_" + os.path.basename(funcfile))

    return undist_file
    
