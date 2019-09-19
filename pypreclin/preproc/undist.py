##########################################################################
# NSAp - Copyright (C) CEA, 2013 - 2019
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
B0 inhomogeneties correction.
"""

# Sytem import
import os
import shutil
import subprocess

# Package import
import pypreclin
from pypreclin.utils.export import ungzip_file
from pypreclin.utils.filetools import average_timeserie

# Third party import
from pyconnectome.wrapper import FSLWrapper
from pyconnectome.utils.segtools import bet2
from pyconnectome.utils.preproctools import topup as _topup
from pyconnectome.utils.preproctools import fsl_prepare_fieldmap
from pyconnectome.configuration import environment, concat_environment
from nipype.interfaces.fsl.preprocess import FUGUE


DIR_MAP = {
    "i": "x",
    "i-": "x-",
    "j": "y",
    "j-": "y-",
    "k": "z",
    "k-": "z-"
}


def mc_undist(funcfile, outdir, spmpath, fslpath):
    """ Slice by slice motion correction: deprecated
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


def topup(blip_up_file, blip_down_file, blip_up_phase_enc_dir,
          blip_down_phase_enc_dir, outdir, fsl_sh, apply_to=None,
          unwarp_direction=None, dwell_time=None, verbose=0):
    """ Bias field correction using topup.

    Note: we do not set the total readout times for both acquisitons, and we
    expect that the readout time is identical for all acquisitions. This 
    parameter simply serves to ensure that the estimated field is correctly
    scaled in Hz, but it doesn't affect the result in terms of undistorted
    images. We set it to 1, meaning that the field will be scaled as 'voxels 
    displacement' instead.

    Parameters
    ----------
    blip_up_file:  str
        path to EPI file acquired in opposite phase enc. directions: blip up.
    blip_down_file:  str
        path to EPI file acquired in opposite phase enc. directions: blip down.
    blip_up_phase_enc_dir: str
        the phase enc. direction for the blip up.
    blip_down_phase_enc_dir: str
        the phase enc. direction for the blip down.
    phase_enc_dirs: list of str
        the phase enc. directions.
    outdir: str
        directory for outputs.
    fsl_sh: str
        path to fsl setup sh file.
    apply_to: str, default None
        apply the topup correction to the EPI volume.
    unwarp_direction: str, default None
        apply topup correction in this direction.
    dwell_time: float, default None
        set the EPI dwell time per phase-encode line, same as echo spacing (ms).
    verbose: int, default 0
        control the verbosity level.

    Returns
    -------
    fieldmap_hz_file: str
        the fieldmap in Hz
    unwarped_epi_file: str
        unwarped EPI file.
    """
    # Call topup
    fieldmap_hz_file, _, _ = _topup(
        b0s=[blip_up_file, blip_down_file],
        phase_enc_dirs=[blip_up_phase_enc_dir, blip_down_phase_enc_dir],
        readout_time=1,
        outroot=outdir,
        apply_to=None,
        fsl_sh=fsl_sh)

    # Convert the fieldmap in rad/s to Hz
    fieldmap_file = fieldmap_hz_file.replace(".nii.gz", "_rads.nii.gz")
    cmd = ["fslmaths", fieldmap_hz_file, "-mul", "6.28", fieldmap_file]
    fslprocess = FSLWrapper(shfile=fsl_sh)
    fslprocess(cmd=cmd)

    # Update env
    fslenv = environment(fsl_sh)
    if (fslenv["FSLDIR"] != os.environ.get("FSLDIR", "")):
        os.environ = concat_environment(os.environ, fslenv)

    # Apply topup
    # Unwarping an input image (fieldmap in rad/s is known)
    unwarp_direction = DIR_MAP[unwarp_direction]
    unwarped_epi_file = os.path.join(outdir, "unwarped_epi_file.nii.gz")
    fugue = FUGUE()
    fugue.inputs.in_file = apply_to
    fugue.inputs.fmap_in_file = fieldmap_file
    fugue.inputs.unwarp_direction = unwarp_direction
    fugue.inputs.unwarped_file = unwarped_epi_file
    fugue.inputs.dwell_time = dwell_time * 1e-3
    fugue.inputs.output_type = "NIFTI_GZ"
    if verbose > 0:
        print(fugue.cmdline)
    returncode = fugue.run()

    return fieldmap_hz_file, unwarped_epi_file


def fugue(epi_file, phase_file, magnitude_file, delta_te, dwell_time,
          unwarp_direction, manufacturer, outdir, fsl_sh, verbose=0):
    """ Unwarping of an EPI image based on fieldmap data using fugue.

    Note: Brain extraction of the EPI image is very important and must be
    tight - that is, it must exclude all non-brain voxels and any voxels with
    only a small partial volume contribution. The reason for this is that
    these areas are normally very noisy in the phase image. The exclusion
    of brain voxels is actually fine and will have no repercussions, since the
    fieldmap is extrapolated beyond this mask.

    Note: If parallel acceleration is used in the EPI acquisition then the
    *effective* echo spacing (dwell time) is the actual echo spacing between
    acquired lines in k-space divided by the acceleration factor.

    Parameters
    ----------
    epi_file: str
        the EPI file to unwarp.
    phase_file: str
        the phase image in the EPI space.
    magnitude_file: str
        the magnitude image in the EPI space.
    delta_te: float
        the echo time difference of the fieldmap sequence - find this out form
        the operator (defaults are *usually* 2.46ms on SIEMENS).
    dwell_time: float
        set the EPI dwell time per phase-encode line, same as echo spacing (ms).
    unwarp_direction: str
        specifies direction of warping: ('x' or 'y' or 'z' or 'x-' or 'y-'
        or 'z-')
    manufacturer: str
        must be SIEMENS.
    outdir: str
        directory for outputs.
    fsl_sh: str
        path to fsl setup sh file.
    verbose: int, default 0
        control the verbosity level.

    Returns
    -------
    magnitude_brain_mask_file: str
        the brain mask.
    vsm_file: str
        voxel shift map file.
    unwarped_epi_file: str
        unwarped EPI file.
    """
    # Extract brain
    avgmagnitude_file = os.path.join(outdir, "avgmag.nii.gz")
    average_timeserie(magnitude_file, avgmagnitude_file, fslconfig=fsl_sh)
    outputs = bet2(avgmagnitude_file, outdir, mask=True, f=0.35, shfile=fsl_sh)
    #outputs = bet2(outputs[0], outdir, mask=True, f=0.6, shfile=fsl_sh)
    magnitude_brain_file, magnitude_brain_mask_file = outputs[:2]

    # Update env
    fslenv = environment(fsl_sh)
    if (fslenv["FSLDIR"] != os.environ.get("FSLDIR", "")):
        os.environ = concat_environment(os.environ, fslenv)

    # Prepare a rad/s fieldmap
    fieldmap_file = os.path.join(outdir, "fieldmap.nii.gz")
    fieldmap_file, fieldmap_hz_file = fsl_prepare_fieldmap(
        manufacturer=manufacturer.upper(),
        phase_file=phase_file,
        brain_magnitude_file=magnitude_brain_file,
        output_file=fieldmap_file,
        delta_te=str(delta_te),
        fsl_sh=fsl_sh)

    # Computing the vsm (unwrapped phase map is known)
    unwarp_direction = DIR_MAP[unwarp_direction]
    vsm_file = os.path.join(outdir, "vsm.nii.gz")
    fugue = FUGUE()
    fugue.inputs.fmap_in_file = fieldmap_file
    fugue.inputs.mask_file = magnitude_brain_mask_file
    fugue.inputs.dwell_time = dwell_time * 1e-3
    fugue.inputs.asym_se_time = delta_te * 1e-3
    fugue.inputs.unwarp_direction = unwarp_direction
    fugue.inputs.save_unmasked_shift = True
    fugue.inputs.shift_out_file = vsm_file
    fugue.inputs.smooth2d = 2
    fugue.inputs.output_type = "NIFTI_GZ"
    fugue.environ = {}
    if verbose > 0:
        print(fugue.cmdline)
    returncode = fugue.run()
    #outputs = returncode.outputs.get()
    #vsm_file = outputs["shift_out_file"]

    # Unwarping an input image (shift map is known)
    unwarped_epi_file = os.path.join(outdir, "unwarped_epi_file.nii.gz")
    fugue = FUGUE()
    fugue.inputs.in_file = epi_file
    fugue.inputs.mask_file = magnitude_brain_mask_file
    fugue.inputs.shift_in_file = vsm_file
    fugue.inputs.unwarp_direction = unwarp_direction
    fugue.inputs.unwarped_file = unwarped_epi_file
    fugue.inputs.output_type = "NIFTI_GZ"
    if verbose > 0:
        print(fugue.cmdline)
    returncode = fugue.run()
    #outputs = returncode.outputs.get()
    #unwarped_epi_file = outputs["unwarped_file"]

    return magnitude_brain_mask_file, vsm_file, unwarped_epi_file


        
