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
import nibabel
import numpy
import shutil

# Hopla import
from hopla.converter import hopla

# Package import
from pypreclin.utils.export import ungzip_file
from pypreclin.utils.export import gzip_file
from pypreclin import preproc

# Pyconnectome
import pyconnectome
from pyconnectome.utils.filetools import apply_mask
from pyconnectome import DEFAULT_FSL_PATH


def timeserie_to_reference(tfile, rindex, outdir, njobs=1, clean_tmp=True):
    """ Register all the fMRI volumes to a reference volume identified by his
    index in the timeserie.

    The registration used is a non-linear regisration.

    Parameters
    ----------
    tfile: str
        the file that containes the timeserie.
    rindex: int
        the reference volume index in the timeserie.
    outdir: str
        the destination folder.
    njobs: int, default 1
        the desired number of parallel job during the registration.
    clean_tmp: bool, default True
        if set, clean the temporary results.

    Returns
    -------
    rfile: str
        the registration result.
    """
    # Check input index
    im = nibabel.load(tfile)
    array = im.get_data()
    if array.ndim != 4:
        raise ValueError("A timeserie (4d volume) is expected.")
    if rindex < 0 or rindex > array.shape[3]:
        raise ValueError(
            "Index '{0}' is out of bound considering the last dimension as "
            "the time dimension '{1}'.".format(rindex, array.shape))

    # Split the timeserie
    tmpdir = os.path.join(outdir, "tmp")
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    moving_images = []
    outdirs = []
    reference_image = None
    for i in range(array.shape[3]):
        _im = nibabel.Nifti1Image(array[..., i], im.affine)
        _outfile = os.path.join(tmpdir, str(i).zfill(4) + ".nii.gz")
        outdirs.append(os.path.join(tmpdir, str(i).zfill(4)))
        if not os.path.isdir(outdirs[-1]):
            os.mkdir(outdirs[-1])
        nibabel.save(_im, _outfile)
        moving_images.append(_outfile)
        if i == rindex:
            reference_image = _outfile

    # Start volume to volume non rigid registration
    scriptdir = os.path.join(os.path.dirname(pyconnectome.__file__), "scripts")
    logfile = os.path.join(tmpdir, "log")
    status, exitcodes = hopla(
        os.path.join(scriptdir, "pyconnectome_ants_register"),
        b="/usr/lib/ants",
        o=outdirs,
        i=moving_images,
        r=reference_image,
        w=1,
        D=3,
        G=0.2,
        J=1,
        N=True,
        B=True,
        v=2,
        hopla_iterative_kwargs=["o", "i"],
        hopla_cpus=njobs,
        hopla_logfile=logfile,
        hopla_verbose=1)
    if not (numpy.asarray(exitcodes.values()) == 0).all():
        raise ValueError("The registration failed, check the log "
                         "'{0}'.".format(logfile))

    # Start timeserie concatenation
    timeserie = []
    affine = None
    for path in outdirs:
        _im = nibabel.load(os.path.join(path,  "ants_2WarpToTemplate.nii.gz"))
        if affine is None:
            affine = _im.affine
        elif not numpy.allclose(affine, im.affine):
            raise ValueError("Affine matrices must be the same.")
        data = _im.get_data()
        data.shape += (1, )
        timeserie.append(data)
    registered_array = numpy.concatenate(timeserie, axis=3)
    _im = nibabel.Nifti1Image(registered_array, affine)
    rfile = os.path.join(outdir, "ants_WarpToTemplate.nii.gz")
    nibabel.save(_im, rfile)

    # Clean temporary files if requested
    if clean_tmp:
        shutil.rmtree(tmpdir, ignore_errors=True)

    return rfile


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
    batch = os.path.abspath(batch)

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


    
