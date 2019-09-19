|Python35|_

.. |Python35| image:: https://img.shields.io/badge/python-3.5-blue.svg
.. _Python35: https://badge.fury.io/py/pypreclin



Description
===========

pypreclin is a Python project that provides a collection of Python scripts for
preprocessing MRI preclinical datasets.
This work is made available by a community of people, amoung which the
CEA Neurospin UNATI and CEA NeuroSpin UNICOG laboratories, in particular A. Grigis,
J. Tasserie, and B. Jarraya.

Important links
===============

- Official source code repo: https://github.com/neurospin/pypreclin
- HTML documentation (last stable release): http://neurospin.github.io/pypreclin

Dependencies
============

The required Python dependencies to use the software are:

* numpy
* scipy
* hopla
* pyconnectome
* pyconnectomist
* nipype
* matplotlib
* nibabel
* joblib
* transforms3d
* filelock
* python-pypipe (optional)

This package requires also other softwares:

* FSL: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation
* ANTS: http://stnava.github.io/ANTs
* JIP-align: http://www.nmr.mgh.harvard.edu/~jbm/jip
* SPM12-standalone: http://www.fil.ion.ucl.ac.uk/spm/software/spm12

Install
=======

Make sure you have installed all the dependencies listed above properly.
Further instructions are available at http://neurospin.github.io/pypreclin/generated/installation.html

Using Singularity
=================

Singularity (https://singularity.lbl.gov/) is convenient to deploy complete
pypreclin installations (including all dependencies). Once deployed, it is an
isolated environment with separate filesystem and namespaces for processes,
etc. Please refer to the Singularity documentation for more details.

We provide a Singularity image with pypreclin, which can be found at http://biodev.cea.fr/pypreclin/pypreclin-ubuntu.simg.

To use it install first Singularity. On Debian/Ubuntu the package is called
singularity-container. You can check that Singularity is installed by
typing singularity --help in a terminal.


