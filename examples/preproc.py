"""
pypreclin preprocessing
=======================

Credit: A Grigis & J Tasserie

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

import os
from pypreclin.workflow.preproc_fmri import preproc
from pprint import pprint

#############################################################################
# Load a test dataset
# -------------------
#
# Now load a toy single run dataset with 10 volumes (it will be faster).
# We assume here that in your current working directory you have a functional
# data 'func.nii', a structural data 'anat.nii' and a template
# 'mni-resampled_1by1by1.nii'.
#
# Note: pypreclin needs external softwares. If you want to install them from
# your terminal on an Ubuntu computer you can run the following commands:
#::
#
#   sudo apt install fsl-5.0-complete ants -y
#   sudo ln -s /usr/lib/ants/N4BiasFieldCorrection /usr/bin
#   wget https://www.nitrc.org/frs/download.php/7446/jip-Linux-x86_64.tar.gz
#   tar xvzf jip-Linux-x86_64.tar.gz 
#   rm jip-Linux-x86_64.tar.gz
#
# In order to avoid the installation
# of all dependencies, please use the provided container.

rootdir = os.getcwd()
funcfile = os.path.join(rootdir, "func.nii")
anatfile = os.path.join(rootdir, "anat.nii")
template = os.path.join(rootdir, "mni-resampled_1by1by1.nii")
for path in (funcfile, anatfile, template):
    if not os.path.isfile(path):
        raise ValueError("You must specify the input data '{0}'.".format(path))
sid = "test"
tr = 2.4
outdir = os.path.join(rootdir, "outputs")
if not os.path.isdir(outdir):
    os.mkdir(outdir)
jipdir = os.path.join(rootdir, "jip-Linux-x86_64", "bin")
if not os.path.isdir(jipdir):
    raise ValueError("Download JIP binaries before starting.")
fslconf = "/etc/fsl/fsl.sh"
if not os.path.isfile(fslconf):
    raise ValueError("Install FSL before starting.")

#############################################################################
# Run
# ---
#
# Finally run all the prepocessing steps.

outputs = preproc(
        funcfile=funcfile,
        anatfile=anatfile,
        sid=sid,
        outdir=outdir,
        repetitiontime=tr,
        template=template,
        jipdir=jipdir,
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
        delta_te=None,
        dwell_time=None,
        manufacturer=None,
        blip_files=None,
        blip_enc_dirs=None,
        unwarp_direction=None,
        phase_file=None,
        magnitude_file=None,
        anatorient="RIA",
        funcorient="RIA",
        kernel_size=3,
        fslconfig=fslconf,
        normalization_trf=None,
        coregistration_trf=None,
        recon1=False,
        recon2=False,
        auto=True,
        verbose=2)
pprint(outputs)


