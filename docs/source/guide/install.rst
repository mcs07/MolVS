.. _install:

Installation
============

.. sectionauthor:: Matt Swain <m.swain@me.com>

MolVS supports Python versions 2.7 and 3.5+.

There are a variety of ways to download and install MolVS.

Option 2: Use conda (recommended)
---------------------------------

The easiest and recommended way to install is using conda. `Anaconda Python`_ is a self-contained Python environment
that is particularly useful for scientific applications.

If you don't already have it, start by installing `Miniconda`_, which includes a complete Python distribution and the
conda package manager. Choose the Python 3 version, unless you have a particular reason why you must use Python 2.

To install MolVS, at the command line, run::

    conda config --add channels conda-forge
    conda install molvs


Option 2: Use pip
-----------------

An alternative method is to install using pip::

    pip install molvs

This will download the latest version of MolVS, and place it in your `site-packages` folder so it is automatically
available to all your python scripts.

.. note::

   MolVS requires RDKit, which cannot be installed using pip. On the Mac, you can use Homebrew::

       brew tap mcs07/cheminformatics
       brew install rdkit

   The official RDKit documentation has `installation instructions for a variety of platforms`_.


Option 2: Download the latest release
-------------------------------------

Alternatively, `download the latest release`_ manually and install yourself::

    tar -xzvf MolVS-0.1.1.tar.gz
    cd MolVS-0.1.1
    python setup.py install

The setup.py command will install MolVS in your `site-packages` folder so it is automatically available to all your
python scripts.

Option 3: Clone the repository
------------------------------

The latest development version of MolVS is always `available on GitHub`_. This version is not guaranteed to be
stable, but may include new features that have not yet been released. Simply clone the repository and install as usual::

    git clone https://github.com/mcs07/MolVS.git
    cd MolVS
    python setup.py install

.. _`Anaconda Python`: https://www.continuum.io/anaconda-overview
.. _`Miniconda`: http://conda.pydata.org/miniconda.html
.. _`installation instructions for a variety of platforms`: http://www.rdkit.org/docs/Install.html
.. _`install it using get-pip.py`: http://www.pip-installer.org/en/latest/installing.html
.. _`download the latest release`: https://github.com/mcs07/MolVS/releases/
.. _`available on GitHub`: https://github.com/mcs07/MolVS
