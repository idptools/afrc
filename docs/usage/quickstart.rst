AFRC Quickstart
=========================================================

Installation
************************
To install the Analytical Flory Random Coil (AFRC) package download the package zip (...) and run 

    ``pip install <filename>``

This will install AFRC in a system-wide manner using the version of ``pip`` specified. We recommend using ``pip`` embedded within a ``conda`` environment to avoid any possible dependency issues. If you are unfamiliar with ``conda`` or ``pip`` we highly recommend reading this `this introductory material <http://geohackweek.github.io/Introductory/01-conda-tutorial//>`_.

AFRC requires ``numpy``, although this is dealt with automatically through the ``pip`` installation.

AFRC was developed for Python 3, but should be compatiable with Python 2. 


Usage
************************

AFRC gives you a way to obtain a variety of inter-residue distances from an analytical version of the Flory Random Coil. The input required is an amino acid sequence, and from this distribution and mean values are available through what is essentially an API.

Specifically, the following code outlines the type of information that might be wanted


To do a developmental install, type

``pip install -e .``

Dependencies
************************
