Installation
=========================================================

afrc requires Python 3.10 or later and depends only on NumPy and SciPy.

From PyPI
----------------------------

Install the latest release with::

   pip install afrc

Writing 3D ensembles as XTC trajectories also needs mdtraj, which you can install alongside afrc with::

   pip install "afrc[ensemble]"

Either way, installing afrc also installs the ``afrc-ensemble`` command-line tool (see :doc:`cli`).


From GitHub
----------------------------

Install the latest development version with::

   pip install git+https://github.com/idptools/afrc.git
