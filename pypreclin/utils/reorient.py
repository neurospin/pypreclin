##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Module that can be used to reorient images.
"""

# System import
import os
import numpy
import nibabel


# Global parameters
POSSIBLE_AXES_ORIENTATIONS = [
    "LAI", "LIA", "ALI", "AIL", "ILA", "IAL",
    "LAS", "LSA", "ALS", "ASL", "SLA", "SAL",
    "LPI", "LIP", "PLI", "PIL", "ILP", "IPL",
    "LPS", "LSP", "PLS", "PSL", "SLP", "SPL",
    "RAI", "RIA", "ARI", "AIR", "IRA", "IAR",
    "RAS", "RSA", "ARS", "ASR", "SRA", "SAR",
    "RPI", "RIP", "PRI", "PIR", "IRP", "IPR",
    "RPS", "RSP", "PRS", "PSR", "SRP", "SPR"
]
CORRECTION_MATRIX_COLUMNS = {
    "R": (1, 0, 0),
    "L": (-1, 0, 0),
    "A": (0, 1, 0),
    "P": (0, -1, 0),
    "S": (0, 0, 1),
    "I": (0, 0, -1)
}


def swap_affine(axes):
    """ Build a correction matrix, from the given orientation of axes to RAS.

    Parameters
    ----------
    axes: str (manadtory)
        the given orientation of the axes.

    Returns
    -------
    rotation: array (4, 4)
        the correction matrix.
    """
    rotation = numpy.eye(4)
    rotation[:3, 0] = CORRECTION_MATRIX_COLUMNS[axes[0]]
    rotation[:3, 1] = CORRECTION_MATRIX_COLUMNS[axes[1]]
    rotation[:3, 2] = CORRECTION_MATRIX_COLUMNS[axes[2]]
    return rotation


def reorient_image(in_file, axes="RAS", prefix="swap", output_directory=None):
    """ Rectify the orientation of an image in order to be in the 'RAS'
    coordinate system.

    Parameters
    ----------
    in_file: str (mandatory)
        the input image.
    axes: str (optional, default 'RAS')
        orientation of the original axes X, Y, and Z
        specified with the following convention: L=Left, R=Right,
        A=Anterion, P=Posterior, I=Inferior, S=Superior.
    prefix: str (optional, default 'swap')
        prefix of the output image.
    output_directory: str (optional, default None)
        the output directory where the rectified image is saved.
        If None use the same directory as the input image.

    Returns
    -------
    out_file: str
        the rectified image.

    Examples
    --------

    >>> from pyfmri.utils.reorient import reorient_image
    >>> rectified_image = reorient_image('image.nii', 'RAS', 's', None)
    """
    # Check the input image exists on the file system
    if not os.path.isfile(in_file):
        raise ValueError("'{0}' is not a valid filename.".format(in_file))

    # Check that the outdir is valid
    if output_directory is not None:
        if not os.path.isdir(output_directory):
            raise ValueError("'{0}' is not a valid directory.".format(
                output_directory))
    else:
        output_directory = os.path.dirname(in_file)

    # Check that a valid input axes is specified
    if axes not in POSSIBLE_AXES_ORIENTATIONS:
        raise ValueError("Wrong '{0}' coordinate system.".format(axes))

    # Get the transformation to the RAS space
    rotation = swap_affine(axes)
    det = numpy.linalg.det(rotation)
    if det != 1:
        print("Rotation matrix determinant must be one not '{0}'.".format(det))

    # Load the image to rectify
    image = nibabel.load(in_file)

    # Get the input image affine transform
    affine = image.get_affine()

    # Apply the rotation to set the image in the RAS coordiante system
    transformation = numpy.dot(rotation, affine)
    image.set_qform(transformation)
    image.set_sform(transformation)

    # Save the rectified image
    basename = os.path.basename(in_file)
    if basename.endswith(".nii"):
        basename += ".gz"
    elif not basename.endswith(".nii.gz"):
        basename += ".nii.gz"
    out_file = os.path.join(output_directory, prefix + basename)
    nibabel.save(image, out_file)

    return out_file


def switch_radiological_neurological(in_file, out_file):
    """ Radiologists like looking at their images with the patient's left
    on the right of the image. If they are looking at a brain image, it is
    as if they were looking at the brain slice from the point of view of
    the patient's feet. Neurologists like looking at brain images with the
    patient's right on the right of the image. This perspective is as if the
    neurologist is looking at the slice from the top of the patient's head.
    The convention is one of image display.

    Parameters
    ----------
    in_file: str (mandatory)
        the input image in RAS convention.
    out_file: str (mandatory)
        the path to the input image with the left/right direction switched.
    """
    im = nibabel.load(in_file)
    switch_im = nibabel.Nifti1Image(
        im.get_data()[::-1, :, :], affine=im.get_affine(),
        header=im.get_header())
    nibabel.save(switch_im, out_file)


def check_orientation(in_files):
    """ Check if the input images have the same orientation.

    Parameters
    ----------
    in_files: str or list of str (mandatory)
        some images.

    Returns
    -------
    test: bool
        true if the orientations are the same.
    """
    if not isinstance(in_files, list):
        in_files = [in_files]
    orients = []
    for path in in_files:
        im = nibabel.load(path)
        orients.append("".join(nibabel.aff2axcodes(im.affine)))
    return all([elem == orients[0] for elem in orients]), orients


def guess_orientation(in_file):
    """ Return an image orientation.

    Parameters
    ----------
    in_files: str or list of str (mandatory)
        an image.

    Returns
    -------
    axes: str
        the image orientation.
    """
    im = nibabel.load(in_file)
    return "".join(nibabel.aff2axcodes(im.affine))
    
