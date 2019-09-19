"""
pypreclin preprocessing
=======================

Credit: A Grigis & J Tasserie

WORK IN PROGRESS

The preprocessing and analysis of nonhuman primate (NHP) magnetic resonance
imaging (MRI) data presents some unique challenges.

1. non standard orientation, sphinx position.
2. strong intensity bias.
3. non brain tissues.

Over the years, we have created our own custom solution that solves these
problems.

This example shows how to preprocess a functional timeserie without
fielmap or reverse phase encoded images. One can inpsect the called function to
apply only a single step.

First checks
------------

In order to test if the 'pypreclin' package is installed on your machine, you
can try to import it.
"""

from pypreclin.workflow.preproc_fmri import prepoc
from pprint import pprint

#############################################################################
# Load a test dataset
# -------------------
#
# Now load a toy single run dataset with 10 volumes.

rootdir = "/tmp"
funcfile = os.path.join(rootdir, "")
anatfile = ""
sid = ""
tr = ""
outdir = ""
template = ""
jipdir = ""

#############################################################################
# Run
# ---
#
# Finally run all the prepocessing steps. In order to avoid the installation
# of all dependencies, please use the provided container.

outputs = preproc(
        funcfile=funcfile,
        anatfile=anatfile,
        sid=sid,
        outdir=outdir,
        repetitiontime=tr,
        template=template,
        jipdir="/jip/bin",
        erase=False,
        resample=False,
        interleaved=False,
        sliceorder="ascending",
        realign_dof=6,
        realign_to_vol=True,
        warp=False,
        warp_njobs=1,
        warp_index=8,
        warp_file=None,
        warp_restrict=[0, 1, 0],
        anatorient="RAS",
        funcorient="RAS",
        kernel_size=3,
        fslconfig="/etc/fsl/5.0.11/fsl.sh",
        normalization_trf=None,
        coregistration_trf=None,
        recon1=False,
        recon2=False,
        auto=True,
        verbose=2)
pprint(outputs)


