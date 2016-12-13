##########################################################################
# NSAp - Copyright (C) CEA, 2013 - 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

# System import
import os
import subprocess

# Module import
from pypreclin.utils.export import ungzip_file
from pypreclin.utils.export import gzip_file
from pypreclin import preproc
from pyconnectome.utils.filetools import apply_mask
from pyconnectome import DEFAULT_FSL_PATH


def jip_align(source_file, target_file, outdir, jipdir, prefix="w",
              fslconfig=DEFAULT_FSL_PATH):
    """ Register a source image to a taget image using the 'jip_align'
    command.

    Parameters
    ----------
    source_file: str (mandatory)
        the source Nifti image.
    target_file: str (mandatory)
        the target Nifti masked image.
    outdir: str (mandatory)
        the destination folder.
    jipdir: str (mandatory)
        the jip binary path.
    prefix: str (optional, default 'w')
        prefix the generated file with this character.
    fslconfig: str (optional)
        the FSL .sh configuration file.

    Returns
    -------
    register_file: str
        the registered image.
    register_mask_file: str
        the registered and masked image.
    native_masked_file: str
        the masked image in the native space.
    """
    # Change current working directory
    cwd = os.getcwd()
    os.chdir(outdir)

    # Only support uncompressed Nifti
    source_file = ungzip_file(source_file, prefix="", outdir=outdir)
    target_file = ungzip_file(target_file, prefix="", outdir=outdir)

    # Create jip environment
    jip_envriron = os.environ
    if "PATH" in jip_envriron:
        jip_envriron["PATH"] = jip_envriron["PATH"] + ":" + jipdir
    else:
        jip_envriron["PATH"] = jipdir

    # Copy source file
    align_file = os.path.join(outdir, "align.com")
    if not os.path.isfile(align_file):
        cmd = ["align", source_file, "-t", target_file]
        print cmd
        subprocess.call(cmd, env=jip_envriron)
        if not os.path.isfile(align_file):
            raise ValueError(
                "No 'align.com' file in '{0}' folder.".format(outdir))

    # Get the apply nonlinear deformation jip batches
    aplly_nonlin_batch = os.path.join(
        os.path.dirname(preproc.__file__), "resources", "apply_nonlin.com")
    aplly_inv_nonlin_batch = os.path.join(
        os.path.dirname(preproc.__file__), "resources", "apply_inv_nonlin.com")

    # Resample the source image
    register_file = os.path.join(outdir, prefix + os.path.basename(source_file))
    cmd = ["jip", aplly_nonlin_batch, source_file, register_file]
    subprocess.call(cmd, env=jip_envriron)

    # Apply mask
    if os.path.isfile(register_file + ".gz"):
        os.remove(register_file + ".gz")
    register_mask_fileroot = os.path.join(
        outdir, "m" + prefix + os.path.basename(source_file).split(".")[0])
    register_mask_file = apply_mask(
        input_file=register_file,
        output_fileroot=register_mask_fileroot,
        mask_file=target_file,
        fslconfig=fslconfig)

    # Send back masked image to original space
    register_mask_file = ungzip_file(register_mask_file, prefix="",
                                     outdir=outdir)
    native_masked_file = os.path.join(
        outdir, "n" + prefix + os.path.basename(register_mask_file))
    cmd = ["jip", aplly_inv_nonlin_batch, register_mask_file,
           native_masked_file]
    subprocess.call(cmd, env=jip_envriron)              

    # Restore current working directory and gzip output
    os.chdir(cwd)
    register_file = gzip_file(
        register_file, prefix="", outdir=outdir,
        remove_original_file=True)
    register_mask_file = gzip_file(
        register_mask_file, prefix="", outdir=outdir,
        remove_original_file=True)
    native_masked_file = gzip_file(
        native_masked_file, prefix="", outdir=outdir,
        remove_original_file=True)

    return register_file, register_mask_file, native_masked_file, align_file


def apply_jip_align(apply_to_files, align_with, outdir, jipdir, prefix="w",
                    apply_inv=False):
    """ Apply a jip deformation.

    Parameters
    ----------
    apply_to_files: list of str (mandatory)
        a list of image path to apply the computed deformation.
    align_with: str ot list of str (mandatory)
        the alignement file containind the deformation parameters. Expect
        an 'align.com' file.
    outdir: str (mandatory)
        the destination folder.
    jipdir: str (mandatory)
        the jip binary path.
    prefix: str (optional, default 'w')
        prefix the generated file with this character.
    apply_inv: bool (optional, default False)
        if set apply the inverse deformation.

    Returns
    -------
    deformed_files: list of str
        the deformed input files.
    """
    # Check input parameters
    if not isinstance(align_with, list):
        align_with = [align_with]
    if not isinstance(apply_to_files, list):
        raise ValueError("The 'apply_to_files' function parameter must "
                         "contains a list of files.")
    for path in apply_to_files + align_with:
        if not os.path.isfile(path):
            raise ValueError("'{0}' is not a valid file.".format(path))

    # Get the apply nonlinear deformation jip batches
    if apply_inv:
        batch = os.path.join(os.path.dirname(preproc.__file__), "resources",
                             "apply_inv_nonlin.com")
    else:
        batch = os.path.join(os.path.dirname(preproc.__file__), "resources",
                             "apply_nonlin.com")

    # Create jip environment
    jip_envriron = os.environ
    if "PATH" in jip_envriron:
        jip_envriron["PATH"] = jip_envriron["PATH"] + ":" + jipdir
    else:
        jip_envriron["PATH"] = jipdir

    # Apply the jip deformation
    deformed_files = []
    cwd = os.getcwd()
    for path in apply_to_files:
        extra_file = ungzip_file(path, prefix="", outdir=outdir)
        deformed_file = os.path.join(
            outdir, prefix + os.path.basename(extra_file))
        for align_file in align_with:
            os.chdir(os.path.dirname(align_file))
            cmd = ["jip", batch, extra_file, deformed_file]
            subprocess.call(cmd, env=jip_envriron)
            extra_file = deformed_file
        deformed_file = gzip_file(
            deformed_file, prefix="", outdir=outdir, remove_original_file=True)
        deformed_files.append(deformed_file)
    os.chdir(cwd)

    return deformed_files


    
