|Python27|_

.. |Python27| image:: https://img.shields.io/badge/python-2.7-blue.svg
.. _Python27: https://badge.fury.io/py/python-pySAP



Description
===========

pyPreClin is a Python project that provides a collection of python scripts for
preprocessing MRI preclinical datasets.
This work is made available by a community of people, amoung which the
CEA Neurospin UNATI and CEA NeuroSpin UNICOG laboratories, in particular A. Grigis,
J. Tasserie, and B. Jarraya.

Important links
===============

- Official source code repo: xxx
- HTML documentation (last stable release): xxx

Dependencies
============

The required Python dependencies to use the software are:

* numpy
* hopla
* bredala
* pyconnectome
* pyconnectomist
* nipype
* matplotlib
* pypreprocess

This package requires also other softwares:

* FSL: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation
* ANTS: http://stnava.github.io/ANTs
* JIP-align: http://www.nmr.mgh.harvard.edu/~jbm/jip
* SPM12-standalone: http://www.fil.ion.ucl.ac.uk/spm/software/spm12

Install
=======

Make sure you have installed all the dependencies listed above properly.
Further instructions are available at xxx

Using Singularity
=================

Singularity (https://singularity.lbl.gov/) is a convenient to deploy complete
pyPreClin installations (including all dependencies). Once deployed, it is an
isolated environment with separate filesystem and namespaces for processes,
etc. Please refer to the Singularity documentation at for more details.

We provide a Singularity image with pyPreClin, which can be found at xxx.

To use it install first Singularity. On Debian/Ubuntu the package is called
singularity-container. You can check that Singularity is installed by
typing singularity --help in a terminal.


