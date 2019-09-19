##########################################################################
# NSAp - Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################


"""
Tools to play with affine matrices.
"""


# Package import
from transforms3d.quaternions import mat2quat, quat2axangle

# Third party import
import numpy
import scipy.linalg

# Global parameters
MAX_ANGLE = 1e10 * 2 * numpy.pi
SMALL_ANGLE = 1e-30


def matrix44_to_vec12(affine):
    """ Decompose a 4x4 matrix describing a homogenous affine transform
    into a 12-sized vector composed of the translations, rotations, zooms,
    shears.

    Decomposing a matrix into simple transformations* by Spencer
    W. Thomas, pp 320-323 in *Graphics Gems II*, James Arvo (editor),
    Academic Press, 1991, ISBN: 0120644819.

    Parameters
    ----------
    affine: array (4, 4)
        the input affine matrix.

    Returns
    -------
    vec12: array (12, )
        the 12-sized vector of natural affine parameters.
    """
    # Initilize the output vector
    vec12 = numpy.zeros((12, ))

    # Extract translation
    vec12[0:3] = affine[:3, 3]

    # Extract x-scale and norm: for this, take the length of the first column
    # vectors e1
    RZS = affine[:-1, :-1]
    e1, e2, e3 = RZS.T
    sx = numpy.sqrt(numpy.sum(e1**2))
    e1 /= sx

    # Extract shear: orthogonalize e2 with respect to e1
    sx_sxy = numpy.dot(e1, e2)
    e2 -= sx_sxy * e1
    sxy = sx_sxy / sx

    # Extract y-scale and norm
    sy = numpy.sqrt(numpy.sum(e2**2))
    e2 /= sy

    # Extract shear: orthogonalize e3 with respect to e1 and e2
    sx_sxz = numpy.dot(e1, e3)
    sy_syz = numpy.dot(e2, e3)
    e3 -= (sx_sxz * e1 + sy_syz * e2)
    sxz = sx_sxz / sx
    syz = sy_syz / sy

    # Extract z-scale and norm
    sz = numpy.sqrt(numpy.sum(e3**2))
    e3 /= sz

    # Extract rotation
    rotation = numpy.array([e1, e2, e3]).T
    if numpy.linalg.det(rotation) < 0:
        sx *= -1
        rotation[:, 0] *= -1
    
    # Extract rotation angles
    angles = rotation_mat2vec(rotation)

    # Construct output vector
    vec12[3:6] = angles
    vec12[6:9] = numpy.array([sx, sy, sz])
    vec12[9:12] = numpy.array([sxy, sxz, syz])

    return vec12


def rotation_mat2vec(rotation):
    """ Rotation vector from rotation matrix.

    Parameters
    ----------
    rotation: array (3, 3)
        rotation matrix.

    Returns
    -------
    vec: array (3, )
        rotation vector, where norm of vec is the angle theta, and the
        axis of rotation is given by vec / theta.
    """
    ax, angle = quat2axangle(mat2quat(rotation))
    return ax * angle


def rotation_vec2mat(r):
    """ The rotation matrix is given by the Rodrigues formula:

    .. math::

        R = Id + sin(theta)*Sn + (1-cos(theta))*Sn^2

    with:

    .. math::

               0  -nz  ny
        Sn =   nz   0 -nx
              -ny  nx   0

    where n = r / ||r||

    In case the angle ||r|| is very small, the above formula may lead
    to numerical instabilities. We instead use a Taylor expansion
    around theta=0:

    .. math::

        R = I + sin(theta)/tetha Sr + (1-cos(theta))/teta2 Sr^2

    leading to:

    .. math::

        R = I + (1-theta2/6)*Sr + (1/2-theta2/24)*Sr^2

    To avoid numerical instabilities, an upper threshold is applied to
    the angle. It is chosen to be a multiple of 2*pi, hence the
    resulting rotation is then the identity matrix. This strategy warrants
    that the output matrix is a continuous function of the input vector.

    Parameters
    ----------
    vec: array (3, )
        rotation vector, where norm of vec is the angle theta, and the
        axis of rotation is given by vec / theta.

    Returns
    -------
    rotation: array (3, 3)
        rotation matrix.
    """
    theta = numpy.sqrt(numpy.sum(r ** 2))
    if theta > MAX_ANGLE:
        return numpy.eye(3)
    elif theta > SMALL_ANGLE:
        n = r / theta
        Sn = numpy.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
        R = numpy.eye(3) + numpy.sin(theta) * Sn\
            + (1 - numpy.cos(theta)) * numpy.dot(Sn, Sn)
    else:
        Sr = numpy.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])
        theta2 = theta * theta
        R = numpy.eye(3) + (1 - theta2 / 6.) * Sr\
            + (.5 - theta2 / 24.) * numpy.dot(Sr, Sr)
    return R

