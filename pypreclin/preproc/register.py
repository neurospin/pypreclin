##########################################################################
# NSAp - Copyright (C) CEA, 2013 - 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Registration utilities.
"""

# System import
import os
import subprocess
import nibabel
import numpy
import shutil
import multiprocessing

# Hopla import
from hopla.converter import hopla

# Package import
from pypreclin.utils.reorient import check_orientation
from pypreclin.utils.export import ungzip_file
from pypreclin.utils.export import gzip_file
from pypreclin.utils.reorient import swap_affine
from pypreclin import preproc

# Third party import
import pyconnectome
from pyconnectome.utils.filetools import apply_mask
from pyconnectome import DEFAULT_FSL_PATH
from pyconnectome.utils.regtools import flirt
from nipype.interfaces.ants import AffineInitializer
from nipype.interfaces.ants import ApplyTransforms
from transforms3d.affines import decompose44
from transforms3d.euler import mat2euler
from scipy.io import loadmat
import nibabel


def timeserie_to_reference(tfile, outdir, rindex=None,
                           restrict_deformation=(1, 1, 1), rfile=None, njobs=1,
                           clean_tmp=True):
    """ Register all the fMRI volumes to a reference volume identified by his
    index in the timeserie.

    The registration used is a non-linear regisration.

    Parameters
    ----------
    tfile: str
        the file that containes the timeserie.
    outdir: str
        the destination folder.
    rindex: int, default None
        the reference volume index in the timeserie.
    restrict_deformation: 3-uplet
        restrict the deformation in the given axis.
    rfile: str, default None
        the reference volume.
    njobs: int, default 1
        the desired number of parallel job during the registration.
    clean_tmp: bool, default True
        if set, clean the temporary results.

    Returns
    -------
    resfile: str
        the registration result.
    """
    # Check input index
    im = nibabel.load(tfile)
    array = im.get_data()
    if array.ndim != 4:
        raise ValueError("A timeserie (4d volume) is expected.")
    reference_image = None
    if rindex is not None:
        if rindex < 0 or rindex > array.shape[3]:
            raise ValueError(
                "Index '{0}' is out of bound considering the last dimension "
                "as the time dimension '{1}'.".format(rindex, array.shape))
    elif rfile is not None:
        reference_image = rfile
    else:
        raise ValueError("You need to specify a reference file or index.")

    # Split the timeserie
    tmpdir = os.path.join(outdir, "tmp")
    if not os.path.isdir(tmpdir):
        os.mkdir(tmpdir)
    moving_images = []
    outdirs = []
    for i in range(array.shape[3]):
        _im = nibabel.Nifti1Image(array[..., i], im.affine)
        _outfile = os.path.join(tmpdir, str(i).zfill(4) + ".nii.gz")
        outdirs.append(os.path.join(tmpdir, str(i).zfill(4)))
        if not os.path.isdir(outdirs[-1]):
            os.mkdir(outdirs[-1])
        nibabel.save(_im, _outfile)
        moving_images.append(_outfile)
        if reference_image is None and i == rindex:
            reference_image = _outfile

    # Start volume to volume non rigid registration
    scriptdir = os.path.join(os.path.dirname(pyconnectome.__file__), "scripts")
    script = os.path.join(scriptdir, "pyconnectome_ants_register")
    python_cmd = "python3"
    if not os.path.isfile(script):
        script = "pyconnectome_ants_register"
        python_cmd = None
    logfile = os.path.join(tmpdir, "log")
    if os.path.isfile(logfile):
        os.remove(logfile)
    status, exitcodes = hopla(
        script,
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
        R=list(restrict_deformation),
        V=2,
        hopla_iterative_kwargs=["o", "i"],
        hopla_cpus=njobs,
        hopla_logfile=logfile,
        hopla_python_cmd=python_cmd,
        hopla_use_subprocess=True,
        hopla_verbose=1)
    if not (numpy.asarray(list(exitcodes.values())) == 0).all():
        raise ValueError("The registration failed, check the log "
                         "'{0}'.".format(logfile))

    # Start timeserie concatenation
    timeserie = []
    affine = None
    for path in outdirs:
        sid = os.path.basename(path)
        _im = nibabel.load(os.path.join(
            path,  "ants_2WarpToTemplate_{0}.nii.gz".format(sid)))
        if affine is None:
            affine = _im.affine
        elif not numpy.allclose(affine, im.affine, atol=1e-3):
            raise ValueError("Affine matrices must be the same: {0} - "
                             "{1}.".format(outdirs[0], path))
        data = _im.get_data()
        data.shape += (1, )
        timeserie.append(data)
    registered_array = numpy.concatenate(timeserie, axis=3)
    _im = nibabel.Nifti1Image(registered_array, affine)
    resfile = os.path.join(outdir, "ants_WarpToTemplate.nii.gz")
    nibabel.save(_im, resfile)

    # Clean temporary files if requested
    if clean_tmp:
        shutil.rmtree(tmpdir, ignore_errors=True)

    return resfile


def jip_align(source_file, target_file, outdir, jipdir, prefix="w",
              auto=False, non_linear=False, fslconfig=DEFAULT_FSL_PATH):
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
    auto: bool (optional, default False)
        if set control the JIP window with the script.
    non_linear: bool (optional, default False)
        in the automatic mode, decide or not to compute the non-linear
        deformation field.
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
    # Check input image orientation: must be the same
    same_orient, orients = check_orientation([source_file, target_file])
    if not same_orient:
        print(
            "[WARNING] Source file '{0}' ({2}) and taget file '{1}' ({3}) "
            "must have the same orientation for JIP to work properly.".format(
                source_file, target_file, orients[0], orients[1]))

    # Try to init the affine deformation
    # ToDo: find a way to init the JIP deformation
    if False:
        nb_cpus = multiprocessing.cpu_count()
        init_file = os.path.join(outdir, "init_affine.mat")
        init = AffineInitializer()
        init.inputs.fixed_image = target_file
        init.inputs.moving_image = source_file
        init.inputs.num_threads = nb_cpus
        init.inputs.out_file = init_file
        init.run()
        init_source_file = os.path.join(outdir, "init_source.nii.gz")
        at = ApplyTransforms()
        at.inputs.dimension = 3
        at.inputs.input_image = source_file
        at.inputs.reference_image = target_file
        at.inputs.output_image = init_source_file
        at.inputs.transforms = init_file
        at.inputs.interpolation = "BSpline"
        at.inputs.num_threads = nb_cpus
        at.run()
        # From https://github.com/ANTsX/ANTs/wiki/ITK-affine-transform-conversion
        # From https://github.com/netstim/leaddbs/blob/master/helpers/ea_antsmat2mat.m
        init_mat = loadmat(init_file)
        afftransform = init_mat["AffineTransform_double_3_3"]
        m_center = init_mat["fixed"]
        mat = numpy.eye(4)
        mat[:3, :3] = afftransform[:9].reshape(3, 3).T
        m_translation = afftransform[9:]
        m_offset = m_translation + m_center - numpy.dot(mat[:3, :3], m_center)
        mat[:3, 3] = m_offset.squeeze()
        mat = numpy.linalg.inv(mat)
        #rotation = swap_affine(axes="LPS")
        #mat = numpy.dot(rotation, mat)
        ras_to_lps = numpy.ones((4, 4))
        ras_to_lps[2, :2] = -1
        ras_to_lps[:2, 2:] = -1
        mat = mat * ras_to_lps
        print(mat)
        trans, rot, zoom, shear = decompose44(mat)
        euler = numpy.asarray(mat2euler(rot, axes='sxyz')) * 180. / numpy.pi
        print(trans, euler, zoom, shear)
        align_file = os.path.join(outdir, "align.com")
        im = nibabel.load(source_file)
        orig = im.affine[:3, 3]
        print(orig)
        shift = orig - trans
        with open(align_file, "wt") as of:
            of.write("set registration-translations {0} \n".format(
                " ".join([str(elem) for elem in shift])))

    # Change current working directory
    cwd = os.getcwd()
    os.chdir(outdir)

    # Only support uncompressed Nifti
    source_file = ungzip_file(source_file, prefix="", outdir=outdir)
    target_file = ungzip_file(target_file, prefix="", outdir=outdir)

    # Create jip environment
    jip_envriron = os.environ
    jip_envriron["JIP_HOME"] = os.path.dirname(jipdir)
    if "PATH" in jip_envriron:
        jip_envriron["PATH"] = jip_envriron["PATH"] + ":" + jipdir
    else:
        jip_envriron["PATH"] = jipdir

    # Copy source file
    align_file = os.path.join(outdir, "align.com")
    cmd = ["align", source_file, "-t", target_file]
    if auto:
        if non_linear:
            auto_cmd = cmd + ["-L", "111111111111", "-W", "111", "-a"]
        else:
            auto_cmd = cmd + ["-L", "111111000000", "-W", "000", "-A"]
        if os.path.isfile(align_file):
            auto_cmd += ["-I"]
        subprocess.call(auto_cmd, env=jip_envriron)  
    else:
        print(" ".join(cmd))
        subprocess.call(cmd, env=jip_envriron)
    if not os.path.isfile(align_file):
        raise ValueError(
            "No 'align.com' file in '{0}' folder. JIP has probably failed: "
            "'{1}'".format(outdir, " ".join(cmd)))

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
    jip_envriron["JIP_HOME"] = os.path.dirname(jipdir)
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


def check_jip_install(jipdir):
    """ Simple function to test if the JIP software bibaries are present.

    Parameters
    ----------
    jip_dir: str
        the JIP binaries directory.

    Returns
    -------
    status: bool
        the JIP installation status.
    """
    status = True
    for name in ("jip", "align"):
        binary_path = os.path.join(jipdir, name)
        if not os.path.isfile(binary_path):
            status = False
            break
    return status


def resample_image(source_file, target_file, out_file,
                   fslconfig=DEFAULT_FSL_PATH):
    """ Resample a source image to fit witha target image.

    Parameters
    ----------
    source_file: str
        the source Nifti image.
    target_file: str
        the target Nifti image.
    out_file: str
        the destination file.
    fslconfig: str (optional)
        the FSL .sh configuration file.
    """
    flirt(
        in_file=source_file,
        ref_file=target_file,
        init=None, cost="corratio",
        applyxfm=True,
        out=out_file,
        interp="nearestneighbour",
        shfile=fslconfig)
    return out_file
    
