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

Example: we assume that in your current working directory you have a functional
data 'func.nii', a structural data 'anat.nii' and a template
'mni-resampled_1by1by1.nii'. 
::

	sudo apt install fsl-5.0-complete ants virtualenv -y
	sudo ln -s /usr/lib/ants/N4BiasFieldCorrection /usr/bin
	wget https://www.nitrc.org/frs/download.php/7446/jip-Linux-x86_64.tar.gz
	tar xvzf jip-Linux-x86_64.tar.gz 
	rm jip-Linux-x86_64.tar.gz
	virtualenv -p /usr/bin/python3.5 ./env
	. env/bin/activate
	pip install --no-cache-dir pypreclin
	python env/bin/pypreclin_preproc_fmri -h
	mkdir outputs
	python env/bin/pypreclin_preproc_fmri \
		-f func.nii \
		-a anat.nii \
		-s test \
		-o outputs \
		-r 2.4 \
		-t mni-resampled_1by1by1.nii \
		-j jip-Linux-x86_64/bin \
		-NA RIA \
		-NF RIA \
		-C /etc/fsl/fsl.sh

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

Example: we assume that in the '/volatile/pypreclin' direcctory you have a functional
data 'func.nii', a structural data 'anat.nii' and a template
'mni-resampled_1by1by1.nii'.

Warning: no X11, you need to provide the align.com files.

::

	wget http://biodev.cea.fr/pypreclin/pypreclin-ubuntu.simg
	mkdir /volatile/pypreclin/outputs
	sudo apt install singularity-container
	sudo nano /etc/singularity/singularity.conf
		mount home = no 
	singularity run --cleanenv --bind /volatile/pypreclin pypreclin-ubuntu.simg \
	    -V 2 \
	    -o /volatile/pypreclin/outputs \
	    -s test \
	    -f /volatile/pypreclin/func.nii \
	    -a /volatile/pypreclin/anat.nii \
	    -r 2.40 \
	    -t /volatile/pypreclin/mni-resampled_1by1by1.nii \
	    -NA RIA \
	    -NF RIA \
	    -N normalize/align.com \
	    -M coreg/align.com \
	    -C /etc/fsl/5.0/fsl.sh \
	    -j /i2bm/local/jip-Linux-x86_64/bin




