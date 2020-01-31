# -*- coding: utf-8 -*-
##########################################################################
# NSAp - Copyright (C) CEA, 2019
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
fMRI preprocessings using FSL, SPM, JIP, and ANTS.

Steps
-----

1. Slice Timing correction (with cache).
2. B0 inhomogeneties correction (optional, with cache).
3. Reorient images not in RAS coordinate system and reorient images to match
   the orientation of the standard MNI152 template (with cache).
4. Realign: motion correction - adjust for movement between slices (with
   cache).
5. Normalization: warp images to fit to a standard template brain (with cache).
6. Tissues segmentation and spatial intensity variations correction (with
   cache).
7. Coregistration: overlay structural and functional images - link
   functional scans to anatomical scan (with cache).
8. Wrap functional: resample the functional serie and mask the registered
   serie (with cache).
9. Smooth the functional time serie (with cache).
10. SNAPs: compute some snaps assessing the different processing steps (with
    cache).
11. Reporting: generate a QC reporting (with cache).
"""

# System import
from __future__ import print_function
import os
import re
import sys
import ast
import json
import glob
import shutil
import tempfile
from pprint import pprint
from datetime import datetime

# Module import
import pypreclin
from pypreclin import DEFAULT_FSL_PATH
from pypreclin.utils.reorient import reorient_image
from pypreclin.utils.reorient import guess_orientation
from pypreclin.preproc.register import check_jip_install
from pypreclin.preproc.register import jip_align
from pypreclin.preproc.register import apply_jip_align
from pypreclin.preproc.register import timeserie_to_reference
from pypreclin.preproc.register import resample_image
from pypreclin.preproc.undist import topup
from pypreclin.preproc.undist import fugue
from pypreclin.plotting.check_preprocessing import plot_fsl_motion_parameters

# PyConnectome import
from pyconnectome.configuration import environment, concat_environment
from pyconnectome.utils.filetools import fslreorient2std
from pyconnectome.utils.filetools import apply_mask
from pyconnectome.utils.segtools import fast
from pyconnectome.utils.segtools import bet2
from pyconnectome.utils.regtools import mcflirt
from pyconnectome.utils.regtools import flirt
from pyconnectome.plotting.slicer import triplanar

# PyConnectomist import
from pyconnectomist.utils.pdftools import generate_pdf

# Third party import
import nibabel
import pylab as plt
from hopla.converter import hopla
from joblib import Memory as JoblibMemory
from nipype.interfaces import fsl
from nipype.interfaces import ants
from nipype.caching import Memory as NipypeMemory


# Global parameters
STEPS = {
    "slice_timing": "1-SliceTiming",
    "warp": "2-B0Inhomogeneities",
    "reorient": "3-Reorient",
    "realign": "4-Realign",
    "normalization": "5-Normalization",
    "inhomogeneities": "6-Inhomogeneities",
    "coregistration": "7-Coregistration",
    "wrap": "8-Wrap",
    "smooth": "9-Smooth",
    "snaps": "10-Snaps",
    "report": "11-Report"
}


def preproc(
        funcfile,
        anatfile,
        sid,
        outdir,
        repetitiontime,
        template,
        jipdir,
        erase,
        resample,
        interleaved,
        sliceorder,
        realign_dof,
        realign_to_vol,
        warp,
        warp_njobs,
        warp_index,
        warp_file,
        warp_restrict,
        delta_te,
        dwell_time,
        manufacturer,
        blip_files,
        blip_enc_dirs,
        unwarp_direction,
        phase_file,
        magnitude_file,
        anatorient,
        funcorient,
        kernel_size,
        fslconfig,
        normalization_trf,
        coregistration_trf,
        recon1,
        recon2,
        auto,
        verbose):
    """ fMRI preprocessings using FSL, SPM, JIP, and ANTS.
    """
    # TODO: remove when all controls available in pypipe
    if not isinstance(erase, bool):
        erase = eval(erase)
        resample = eval(resample)
        interleaved = eval(interleaved)
        realign_to_vol = eval(realign_to_vol)
        warp = eval(warp)
        recon1 = eval(recon1)
        recon2 = eval(recon2)
        auto = eval(auto)
        warp_restrict = eval(warp_restrict)
        blip_files = None if blip_files == "" else eval(blip_files)
        blip_enc_dirs = eval(blip_enc_dirs)

    # Read input parameters
    funcfile = os.path.abspath(funcfile)
    anatfile = os.path.abspath(anatfile)
    template = os.path.abspath(template)
    jipdir = os.path.abspath(jipdir)
    realign_to_mean = not realign_to_vol
    subjdir = os.path.join(os.path.abspath(outdir), sid)
    cachedir = os.path.join(subjdir, "cachedir")
    outputs = {}
    if erase and os.path.isdir(subjdir):
        shutil.rmtree(subjdir)
    if not os.path.isdir(cachedir):
        os.makedirs(cachedir)
    nipype_memory = NipypeMemory(cachedir)
    joblib_memory = JoblibMemory(cachedir, verbose=verbose)
    def display_outputs(outputs, verbose, **kwargs):
        """ Simple function to display/store step outputs.
        """
        if verbose > 0:
            print("-" * 50)
            for key, val in kwargs.items():
                print("{0}: {1}".format(key, val))
            print("-" * 50)
        outputs.update(kwargs)

    # Check input parameters
    template_axes = guess_orientation(template)
    if template_axes != "RAS":
        raise ValueError("The template orientation must be 'RAS', '{0}' "
                         "found.".format(template_axes))
    check_jip_install(jipdir)
    if sliceorder not in ("ascending", "descending"):
        raise ValueError("Supported slice order are: ascending & descending.")

    # Slice timing
    fslenv = environment(fslconfig)
    if (fslenv["FSLDIR"] != os.environ.get("FSLDIR", "")):
        os.environ = concat_environment(os.environ, fslenv)
    st_dir = os.path.join(subjdir, STEPS["slice_timing"])
    if not os.path.isdir(st_dir):
        os.mkdir(st_dir)
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
    interface = nipype_memory.cache(fsl.SliceTimer)
    returncode = interface(
        in_file=funcfile,
        interleaved=interleaved,
        slice_direction=3,
        time_repetition=repetitiontime,
        index_dir=False if sliceorder=="ascending" else True,
        out_file=os.path.join(
            st_dir, os.path.basename(funcfile).split(".")[0] + ".nii.gz"))
    st_outputs = returncode.outputs.get()
    slice_time_corrected_file = st_outputs["slice_time_corrected_file"]
    display_outputs(
        outputs, verbose, slice_time_corrected_file=slice_time_corrected_file)

    # B0 inhomogeneities: topup or fugue or None + motion induced
    if warp:
        warp_dir = os.path.join(subjdir, STEPS["warp"])
        if not os.path.isdir(warp_dir):
            os.mkdir(warp_dir)
        if blip_files is not None:
            interface = joblib_memory.cache(topup)
            fieldmap_hz_file, unwarped_epi_file = interface(
                blip_up_file=blip_files[0],
                blip_down_file=blip_files[1],
                blip_up_phase_enc_dir=blip_enc_dirs[0],
                blip_down_phase_enc_dir=blip_enc_dirs[1],
                apply_to=slice_time_corrected_file,
                unwarp_direction=unwarp_direction,
                dwell_time=dwell_time,
                outdir=warp_dir,
                fsl_sh=fslconfig)
        elif phase_file is not None:
            interface = joblib_memory.cache(fugue)
            magnitude_brain_mask_file, vsm_file, unwarped_epi_file = interface(
                epi_file=slice_time_corrected_file,
                phase_file=phase_file,
                magnitude_file=magnitude_file,
                delta_te=delta_te,
                dwell_time=dwell_time,
                unwarp_direction=unwarp_direction,
                manufacturer=manufacturer,
                outdir=warp_dir,
                fsl_sh=fslconfig,
                verbose=verbose)
        else:
            unwarped_epi_file = slice_time_corrected_file
        interface = joblib_memory.cache(timeserie_to_reference)
        b0_corrected_file = interface(
            tfile=unwarped_epi_file,
            rindex=warp_index,
            restrict_deformation=warp_restrict,
            rfile=warp_file,
            outdir=warp_dir,
            njobs=warp_njobs,
            clean_tmp=True)
    else:
        b0_corrected_file = slice_time_corrected_file
    display_outputs(
        outputs, verbose, b0_corrected_file=b0_corrected_file)

    # Reorient images not in RAS coordinate system and
    # reorient images to match the orientation of the standard MNI152 template.
    reorient_dir = os.path.join(subjdir, STEPS["reorient"])
    if not os.path.isdir(reorient_dir):
        os.mkdir(reorient_dir)
    reoriented_funcfile = b0_corrected_file
    reoriented_anatfile = anatfile
    interface = joblib_memory.cache(reorient_image)
    if funcorient != "RAS":
        reoriented_funcfile = interface(
            b0_corrected_file,
            axes=funcorient,
            prefix="o",
            output_directory=reorient_dir)
    if anatorient != "RAS":
        reoriented_anatfile = interface(
            anatfile,
            axes=anatorient,
            prefix="o",
            output_directory=reorient_dir)
    standard_funcfile = os.path.join(
        reorient_dir, "d" + os.path.basename(reoriented_funcfile).split(".")[0])
    standard_anatfile = os.path.join(
        reorient_dir, "d" + os.path.basename(reoriented_anatfile).split(".")[0])
    interface = joblib_memory.cache(fslreorient2std)
    standard_funcfile = interface(
        reoriented_funcfile,
        standard_funcfile,
        fslconfig=fslconfig)
    standard_anatfile = interface(
        reoriented_anatfile,
        standard_anatfile,
        fslconfig=fslconfig)
    display_outputs(
        outputs, verbose, standard_funcfile=standard_funcfile,
        standard_anatfile=standard_anatfile)

    # Downsample template
    if resample:
        template = resample_image(
            source_file=template,
            target_file=standard_anatfile,
            out_file=os.path.join(subjdir, "template.nii.gz"),
            fslconfig=fslconfig)

    # Realign
    realign_dir = os.path.join(subjdir, STEPS["realign"])
    if not os.path.isdir(realign_dir):
        os.mkdir(realign_dir)
    realign_func_rootfile = os.path.join(
        realign_dir, "r" + os.path.basename(standard_funcfile).split(".")[0])
    interface = joblib_memory.cache(mcflirt)
    realign_funcfile, realign_func_meanfile, realign_func_parfile = interface(
        in_file=standard_funcfile,
        out_fileroot=realign_func_rootfile,
        cost="normcorr",
        bins=256,
        dof=realign_dof,
        refvol=warp_index,
        reffile=warp_file,
        reg_to_mean=realign_to_mean,
        mats=True,
        plots=True,
        verbose=verbose,
        shfile=fslconfig)
    display_outputs(
        outputs, verbose, realign_funcfile=realign_funcfile,
        realign_func_meanfile=realign_func_meanfile)

    # Early stop detected
    if recon1:
        print("[warn] User requested a processing early stop. Remove the 'recon1' "
              "option to resume.")
        return outputs

    # Normalization.
    normalization_dir = os.path.join(subjdir, STEPS["normalization"])
    if not os.path.isdir(normalization_dir):
        os.mkdir(normalization_dir)
    if normalization_trf is not None:
        shutil.copyfile(
            normalization_trf, os.path.join(normalization_dir, "align.com"))
    interface = joblib_memory.cache(jip_align)
    (register_anatfile, register_anat_maskfile,
     native_anat_maskfile, align_normfile) = interface(
        source_file=standard_anatfile,
        target_file=template,
        outdir=normalization_dir,
        jipdir=jipdir,
        prefix="w",
        auto=auto,
        non_linear=True,
        fslconfig=fslconfig)
    display_outputs(
        outputs, verbose, register_anatfile=register_anatfile,
        register_anat_maskfile=register_anat_maskfile,
        native_anat_maskfile=native_anat_maskfile, align_normfile=align_normfile)

    # Tissues segmentation and spatial intensity variations correction.
    inhomogeneities_dir = os.path.join(subjdir, STEPS["inhomogeneities"])
    if not os.path.isdir(inhomogeneities_dir):
        os.mkdir(inhomogeneities_dir)
    biascorrected_anatfile = os.path.join(
        inhomogeneities_dir, "n4_" + os.path.basename(native_anat_maskfile))
    bias_anatfile = os.path.join(
        inhomogeneities_dir, "n4field_" + os.path.basename(native_anat_maskfile))
    n4 = ants.N4BiasFieldCorrection()
    interface = nipype_memory.cache(ants.N4BiasFieldCorrection)
    returncode = interface(
        dimension=3,
        input_image=native_anat_maskfile,
        bspline_fitting_distance=200,
        shrink_factor=2,
        n_iterations=[50,50,40,30],
        convergence_threshold=1e-6,
        output_image=biascorrected_anatfile,
        bias_image=bias_anatfile)
    n4_outputs = returncode.outputs.get()
    display_outputs(
        outputs, verbose, bias_anatfile=bias_anatfile,
        biascorrected_anatfile=biascorrected_anatfile)

    # Coregistration.
    coregistration_dir = os.path.join(subjdir, STEPS["coregistration"])
    if not os.path.isdir(coregistration_dir):
        os.mkdir(coregistration_dir)
    if coregistration_trf is not None:
        shutil.copy(coregistration_trf, coregistration_dir)
    interface = joblib_memory.cache(jip_align)
    (register_func_meanfile, register_func_mean_maskfile,
     native_func_mean_maskfile, align_coregfile) = interface(
        source_file=realign_func_meanfile,
        target_file=biascorrected_anatfile,
        outdir=coregistration_dir,
        jipdir=jipdir,
        prefix="w",
        auto=auto,
        non_linear=False,
        fslconfig=fslconfig)
    display_outputs(
        outputs, verbose, register_func_meanfile=register_func_meanfile,
        register_func_mean_maskfile=register_func_mean_maskfile,
        native_func_mean_maskfile=native_func_mean_maskfile,
        align_coregfile=align_coregfile)

    # Early stop detected
    if recon2:
        print("[warn] User requested a processing early stop. Remove the 'recon2' "
              "option to resume.")
        return outputs

    # Wrap functional: resample the functional serie and mask the registered serie.
    wrap_dir = os.path.join(subjdir, STEPS["wrap"])
    if not os.path.isdir(wrap_dir):
        os.mkdir(wrap_dir)
    interface = joblib_memory.cache(apply_jip_align)
    deformed_files = interface(
        apply_to_files=[realign_funcfile],
        align_with=[align_coregfile, align_normfile],
        outdir=wrap_dir,
        jipdir=jipdir,
        prefix="w",
        apply_inv=False)
    register_funcfile = deformed_files[0]
    register_func_mask_fileroot = os.path.join(
        wrap_dir, "m" + os.path.basename(register_funcfile).split(".")[0])
    interface = joblib_memory.cache(apply_mask)
    register_func_maskfile = interface(
        input_file=register_funcfile,
        output_fileroot=register_func_mask_fileroot,
        mask_file=template,
        fslconfig=fslconfig)
    display_outputs(
        outputs, verbose, register_funcfile=register_funcfile,
        register_func_maskfile=register_func_maskfile)

    # Smooth the functional serie.
    smooth_dir = os.path.join(subjdir, STEPS["smooth"])
    if not os.path.isdir(smooth_dir):
        os.mkdir(smooth_dir)
    interface = nipype_memory.cache(fsl.Smooth)
    returncode = interface(
        in_file=register_func_maskfile,
        fwhm=kernel_size,
        output_type="NIFTI",
        smoothed_file=os.path.join(
            smooth_dir,
            "smooth_" + os.path.basename(register_func_maskfile).split(".")[0] +
            ".nii"))
    smooth_outputs = returncode.outputs.get()
    smoothed_file = smooth_outputs["smoothed_file"]
    display_outputs(outputs, verbose, smoothed_file=smoothed_file)

    # Copy the results to the root directory: use Nifti format.
    nibabel.load(smoothed_file).to_filename(
        os.path.join(subjdir, "sMNI.nii"))
    nibabel.load(register_func_maskfile).to_filename(
        os.path.join(subjdir, "MNI.nii"))
    nibabel.load(register_anat_maskfile).to_filename(
        os.path.join(subjdir, "anat.nii"))

    # Compute some snaps assessing the different processing steps.
    snapdir = os.path.join(subjdir, STEPS["snaps"])
    if not os.path.isdir(snapdir):
        os.mkdir(snapdir)
    interface = joblib_memory.cache(triplanar)
    # > generate coregistration plot
    coregister_fileroot = os.path.join(snapdir, "coregister")
    coregister_file = interface(
        input_file=register_func_meanfile,
        output_fileroot=coregister_fileroot,
        overlays=[standard_anatfile],
        overlays_colors=None,
        contours=True,
        edges=False,
        overlay_opacities=[0.7],
        resolution=300)
    # > generate normalization plot
    normalize_fileroot = os.path.join(snapdir, "normalization")
    normalize_file = interface(
        input_file=register_anatfile,
        output_fileroot=normalize_fileroot,
        overlays=[template],
        overlays_colors=None,
        contours=True,
        edges=False,
        overlay_opacities=[0.7],
        resolution=300)
    # > generate a motion parameter plot
    interface = joblib_memory.cache(plot_fsl_motion_parameters)
    realign_motion_file = os.path.join(snapdir, "realign_motion_parameters.png")
    interface(realign_func_parfile, realign_motion_file)
    display_outputs(
        outputs, verbose, normalize_file=normalize_file,
        realign_motion_file=realign_motion_file, coregister_file=coregister_file)

    # Generate a QC reporting
    reportdir = os.path.join(subjdir, STEPS["report"])
    reportfile = os.path.join(reportdir, "QC_preproc_{0}.pdf".format(sid))
    if not os.path.isdir(reportdir):
        os.mkdir(reportdir)
    interface = joblib_memory.cache(generate_pdf)
    tic = datetime.now()
    date = "{0}-{1}-{2}".format(tic.year, tic.month, tic.day)
    interface(
        datapath=snapdir,
        struct_file=os.path.join(
            os.path.abspath(os.path.dirname(pypreclin.__file__)), "utils",
            "resources", "pypreclin_qcpreproc.json"),
        author="NeuroSpin",
        client="-",
        poweredby="FSL-SPM-Nipype-JIP",
        project="-",
        timepoint="-",
        subject=sid,
        date=date,
        title="fMRI Preprocessing QC Reporting",
        filename=reportfile,
        pagesize=None,
        left_margin=10,
        right_margin=10,
        top_margin=20,
        bottom_margin=20,
        show_boundary=False,
        verbose=0)
    display_outputs(outputs, verbose, reportfile=reportfile)

    return outputs


def preproc_multi(
        bidsdir, template, jipdir, fslconfig, auto, resample, anatorient="RAS",
        funcorient="RAS", njobs=1, simage=None, shome=None, sbinds=None,
        verbose=0):
    """ Perform the FMRI preprocessing on a BIDS organized directory in
    parallel (without FUGUE or TOPUP).
    This function can be called with a singularity image that contains all the
    required software.

    If a 'jip_trf' directory is available in the session directory, the code
    will use the available JIP transformation.

    Parameters
    
    bidsdir: str
        the BIDS organized directory.
    template: str
        the path to the template in RAS coordiante system.
    jipdir: str
        the jip software binary path.
    fslconfig: str
        the FSL configuration script.
    auto: bool
        control the JIP window with the script.
    resample: bool
        if set resample the input template to fit the anatomical image.
    anatorient: str, default "RAS"
        the input anatomical image orientation.
    funcorient: str, default "RAS"
        the input functional image orientation.
    njobs: int, default 1
        the number of parallel jobs.
    simage: simg, default None
        a singularity image.
    shome: str, default None
        a home directory for the singularity image.
    sbinds: str, default None
        shared directories for the singularity image.
    verbose: int, default 0
        control the verbosity level.
    """
    # TODO: remove when all controls available in pypipe
    if not isinstance(auto, bool):
        auto = eval(auto)
        resample = eval(resample)
        sbinds = eval(sbinds)

    # Parse the BIDS directory
    desc_file = os.path.join(bidsdir, "dataset_description.json")
    if not os.path.isfile(desc_file):
        raise ValueError(
            "Expect '{0}' in a BIDS organized directory.".format(desc_file))
    with open(desc_file, "rt") as of:
        desc = json.load(of)
    name = desc["Name"]
    if verbose > 0:
        print("Processing dataset '{0}'...".format(name))
    dataset = {name: {}}
    funcfiles, anatfiles, subjects, outputs, trs, encdirs, normfiles, coregfiles = [
        [] for _ in range(8)]
    outdir = os.path.join(bidsdir, "derivatives", "preproc")
    for sesdir in glob.glob(os.path.join(bidsdir, "sub-*", "ses-*")):
        split = sesdir.split(os.sep)
        sid, ses = split[-2:]
        _anatfiles = glob.glob(os.path.join(sesdir, "anat", "sub-*T1w.nii.gz"))
        _funcfiles = glob.glob(os.path.join(sesdir, "func", "sub-*bold.nii.gz"))
        if len(_anatfiles) != 1:
            if verbose > 0:
                print("Skip session '{0}': no valid anat.".format(sesdir))
            continue
        if len(_funcfiles) == 0:
            if verbose > 0:
                print("Skip session '{0}': no valid func.".format(sesdir))
            continue
        descs = []
        for path in _funcfiles:
            runs = re.findall(".*_(run-[0-9]*)_.*", path)
            if len(runs) != 1:
                if verbose > 0:
                    print("Skip path '{0}': not valid name.".format(path))
                continue
            sesdir = os.path.join(outdir, sid, ses)
            rundir = os.path.join(sesdir, runs[0])
            if not os.path.isdir(rundir):
                os.makedirs(rundir)
            desc_file = path.replace(".nii.gz", ".json")
            if not os.path.isfile(desc_file):
                break
            with open(desc_file, "rt") as of:
                _desc = json.load(of)
            if _desc["PhaseEncodingDirection"].startswith("i"):
                warp_restrict = [1, 0, 0]
            elif _desc["PhaseEncodingDirection"].startswith("j"):
                warp_restrict = [0, 1, 0]
            elif _desc["PhaseEncodingDirection"].startswith("k"):
                warp_restrict = [0, 0, 1]
            else:
                raise ValueError(
                    "Unknown encode phase direction : {0}...".format(
                        _desc["PhaseEncodingDirection"]))
            jip_normalization = os.path.join(
                sesdir, "jip_trf", "Normalization", "align.com")
            if not os.path.isfile(jip_normalization):
                if verbose > 0:
                    print("No JIP normalization align.com file found "
                          "in {0}.".format(sesdir))
                jip_normalization = None
            jip_coregistration = os.path.join(
                sesdir, "jip_trf", "Coregistration", "align.com")
            if not os.path.isfile(jip_coregistration):
                if verbose > 0:
                    print("No JIP coregistration align.com file found "
                          "in {0}.".format(sesdir))
                jip_coregistration = None
            desc = {
                "tr": _desc["RepetitionTime"],
                "warp_restrict": warp_restrict,
                "output": rundir}
            descs.append(desc)
            funcfiles.append(path)
            anatfiles.append(_anatfiles[0])
            subjects.append(sid)
            outputs.append(rundir)
            trs.append(_desc["RepetitionTime"])
            encdirs.append(warp_restrict)
            normfiles.append(jip_normalization)
            coregfiles.append(jip_coregistration)
        if len(_funcfiles) != len(descs):
            if verbose > 0:
                print("Skip session '{0}': no valid desc.".format(sesdir))
            continue
        if sid not in dataset[name]:
            dataset[name][sid] = {}
        dataset[name][sid][ses] = {
            "anat": _anatfiles[0],
            "func": _funcfiles,
            "desc": descs}
    if verbose > 0:
        pprint(dataset)

    # Preprare inputs
    expected_size = len(subjects)
    if expected_size == 0:
        if verbose > 0:
            print("no data to process.")
        return None
    for name, item in (("outputs", outputs), ("funcfiles", funcfiles),
                       ("anatfiles", anatfiles), ("subjects", subjects),
                       ("trs", trs), ("encdirs", encdirs),
                       ("normfiles", normfiles), ("coregfiles", coregfiles)):
        if item is None:
            continue
        if verbose > 0:
            print("[{0}] {1} : {2} ... {3}".format(
                name.upper(), len(item), item[0], item[-1]))
        if len(item) != expected_size:
            raise ValueError("BIDS dataset not aligned.")

    # Run preproc
    scriptdir = os.path.join(os.path.dirname(pypreclin.__file__), "scripts")
    script = os.path.join(scriptdir, "pypreclin_preproc_fmri")
    python_cmd = "python3"
    if not os.path.isfile(script):
        script = "pypreclin_preproc_fmri"
        python_cmd = None
    if simage is not None:
        script = "singularity run --home {0} ".format(
            shome or tempfile.gettempdir())
        sbinds = sbinds or []
        for path in sbinds:
            script += "--bind {0} ".format(path)
        script += simage
        python_cmd = None
    logdir = os.path.join(bidsdir, "derivatives", "logs")
    if not os.path.isdir(logdir):
        os.makedirs(logdir)
    logfile = os.path.join(
        logdir, "pypreclin_{0}.log".format(datetime.now().isoformat()))
    status, exitcodes = hopla(
        script,
        hopla_iterative_kwargs=["f", "a", "s", "o", "r", "WR", "N", "M"],
        f=funcfiles,
        a=anatfiles,
        s=subjects,
        o=outputs,
        r=trs,
        t=template,
        j=jipdir,
        resample=resample,
        W=True,
        WN=1,
        WR=encdirs,
        NA=anatorient,
        NF=funcorient,
        C=fslconfig,
        N=normfiles,
        M=coregfiles,
        A=auto,
        hopla_python_cmd=python_cmd,
        hopla_use_subprocess=True,
        hopla_verbose=1,
        hopla_cpus=njobs,
        hopla_logfile=logfile)

    return {"logfile": logfile}
    
