.. _install:

Installation
============

.. sectionauthor:: Matt Swain <m.swain@me.com>

MolVS supports Python versions 2.7, 3.4 and 3.5.

.. note::

   MolVS requires RDKit to be installed. On Mac OS X, this is easiest with Homebrew::

       brew tap mcs07/cheminformatics
       brew install rdkit

   The official RDKit documentation has `installation instructions for a variety of platforms`_.

There are a variety of ways to download and install MolVS.

Option 1: Use pip (recommended)
-------------------------------

The easiest and recommended way to install is using pip::

    pip install molvs

This will download the latest version of MolVS, and place it in your `site-packages` folder so it is automatically
available to all your python scripts.

.. note::

   If you don't already have pip installed, you can `install it using get-pip.py`_::

          curl -O https://raw.github.com/pypa/pip/master/contrib/get-pip.py
          python get-pip.py

Option 2: Download the latest release
-------------------------------------

Alternatively, `download the latest release`_ manually and install yourself::

    tar -xzvf MolVS-0.1.0.tar.gz
    cd MolVS-0.1.0
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

.. _`installation instructions for a variety of platforms`: http://www.rdkit.org/docs/Install.html
.. _`install it using get-pip.py`: http://www.pip-installer.org/en/latest/installing.html
.. _`download the latest release`: https://github.com/mcs07/MolVS/releases/
.. _`available on GitHub`: https://github.com/mcs07/MolVS
