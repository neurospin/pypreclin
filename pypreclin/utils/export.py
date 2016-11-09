##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
Module that contains export utilities.
"""

# System import
import os
import json
import gzip


def gzip_and_update_outputs(results):
    """ Compress file and update the result outputs accordingly. Note that
    a simplified results structure as produced by the
    'export_interface_results' function is expected.

    Note that the results structure will be updated (ie modified) by the
    function.

    Parameters
    ----------
    results: dict (mandatory)
        a results structure with process name as keys and process results
        as values. Each process values is expected to contain an 'outputs'
        key containing a dictionary mapping of the process outputs.
    """
    for process_name, results in results.items():
        updated_params = []
        if isinstance(results, (list, tuple)):
            for cnt, item in enumerate(results):
               results[cnt] = recursive_gzip(item["outputs"])
            continue
        for param_name, param in results["outputs"].items():
            results["outputs"][param_name] = recursive_gzip(param)


def recursive_gzip(obj):
    """ Recursively find and Gzip files.

    Parameters
    ----------
    obj: object (mandatory)
        a Python object containing files to Gzip.

    Returns
    -------
    gzip_obj: object
        the input object with Gziped files.
    """
    # Stop case: a non iterative structure is reached
    if isinstance(obj, basestring) and os.path.isfile(obj):
        return gzip_file(obj, prefix="", outdir=None,
                         remove_original_file=True)

    # Go deeper
    elif isinstance(obj, (tuple, list)):
        gzip_obj = []
        for item in obj:
            gzip_obj.append(recursive_gzip(obj))
        if isinstance(obj, tuple):
            gzip_obj = tuple(gzip_obj)
        return gzip_obj

    # Default return object
    else:
        return obj


def ungzip_file(fname, prefix="u", outdir=None):
    """ Copy and ungzip the input file.

    Parameters
    ----------
    fname: str (mandatory)
        an input file to ungzip.
    prefix: str (optional, default 'u')
        the prefix of the result file.
    outdir: str (optional, default None)
        the output directory where ungzip file is saved.

    Returns
    -------
    ungzipfname: str
        the returned ungzip file.
    </unit>
    """
    # Check the input file exists on the file system
    if not os.path.isfile(fname):
        raise ValueError("'{0}' is not a valid filename.".format(fname))

    # Check that the outdir is valid
    if outdir is not None:
        if not os.path.isdir(outdir):
            raise ValueError(
                "'{0}' is not a valid directory.".format(outdir))
    else:
        outdir = os.path.dirname(fname)

    # Get the file descriptors
    base, extension = os.path.splitext(fname)
    basename = os.path.basename(base)

    # Ungzip only known extension
    if extension in [".gz"]:

        # Generate the output file name
        basename = prefix + basename
        ungzipfname = os.path.join(outdir, basename)

        # Read the input gzip file
        with gzip.open(fname, "rb") as gzfobj:
            data = gzfobj.read()

        # Write the output ungzip file
        with open(ungzipfname, "wb") as openfile:
            openfile.write(data)

    # Default, unknown compression extension: the input file is returned
    else:
        ungzipfname = fname

    return ungzipfname


def gzip_file(ungzip_file, prefix="g", outdir=None,
              remove_original_file=False):
    """ Gzip an input file and possibly remove the original file.

    Parameters
    ----------
    ungzip_file: str (mandatory)
        an input file to gzip.
    prefix: str (optional, default 'g')
        a prefix that will be concatenated to the produced file basename.
    outdir: str (optional, default None)
        the destination folder where the Gzip file is saved. If this parameter
        is None, the input image folder is considered as an output folder.
    remove_original_file: bool (optiona, default False)
        if True, remove the original file.

    Returns
    -------
    gzip_file: str
        the returned Gzip file.
    """
    # Check the input file exists on the file system
    if not os.path.isfile(ungzip_file):
        raise ValueError("'{0}' is not a valid filename.".format(ungzip_file))

    # Check that the outdir is valid
    if outdir is not None:
        if not os.path.isdir(outdir):
            raise ValueError("'{0}' is not a valid directory.".format(outdir))
    else:
        outdir = os.path.dirname(ungzip_file)

    # Get the file descriptors
    dirname, basename = os.path.split(ungzip_file)
    base, extension = os.path.splitext(basename)

    # Gzip only non compressed file
    if extension not in [".gz", ".mat", ".json", ".txt"]:

        # Generate the output file name
        basename = base + extension + ".gz"
        if prefix:
            basename = prefix + basename
        gzip_file = os.path.join(outdir, basename)

        # Write the output gzip file
        with open(ungzip_file, "rb") as openfile:
            with gzip.open(gzip_file, "w") as gzfobj:
                gzfobj.writelines(openfile)

        # Remove original file if requested
        if remove_original_file:
            os.remove(ungzip_file)

    # Default, the input file is returned
    else:
        gzip_file = ungzip_file

    return gzip_file

