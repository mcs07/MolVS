.. MolVS documentation master file, created by sphinx-quickstart on Thu Apr 24 14:35:38 2014.

MolVS: Molecule Validation and Standardization
==============================================

.. sectionauthor:: Matt Swain <m.swain@me.com>

**MolVS** is a molecule validation and standardization tool, written in Python using the `RDKit chemistry framework`_.

Building a collection of chemical structures from different sources can be difficult due to differing representations,
drawing conventions and mistakes. MolVS can standardize chemical structures to improve data quality, help with
de-duplication and identify relationships between molecules.

There are sensible defaults that make it easy to get started::

    >>> from molvs import standardize_smiles
    >>> standardize_smiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
    '[Na+].O=C([O-])c1ccc(CS(=O)=O)cc1'

Each standardization module is also available separately, allowing the development of custom standardization processes.

Features
--------

- Normalization of functional groups to a consistent format.
- Recombination of separated charges.
- Breaking of bonds to metal atoms.
- Competitive reionization to ensure strongest acids ionize first in partially ionize molecules.
- Tautomer enumeration and canonicalization.
- Neutralization of charges.
- Standardization or removal of stereochemistry information.
- Filtering of salt and solvent fragments.
- Generation of fragment, isotope, charge, tautomer or stereochemistry insensitive parent structures.
- Validations to identify molecules with unusual and potentially troublesome characteristics.

User guide
----------

A step-by-step guide to getting started with MolVS.

.. toctree::
   :maxdepth: 1

   guide/intro
   guide/install
   guide/gettingstarted
   guide/validate
   guide/standardize
   guide/tautomer
   guide/fragment
   guide/charge
   guide/cli
   guide/contributing

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method. This is automatically generated
from the MolVS source code and comments.

.. toctree::
   :maxdepth: 2

   api


Useful links
------------

- `MolVS on GitHub`_
- `MolVS on PyPI`_
- `Issue tracker`_
- `Release history`_
- `MolVS Travis CI`_

.. _`RDKit chemistry framework`: http://www.rdkit.org
.. _`MolVS on GitHub`: https://github.com/mcs07/MolVS
.. _`MolVS on PyPI`: https://pypi.python.org/pypi/MolVS
.. _`Issue tracker`: https://github.com/mcs07/MolVS/issues
.. _`Release history`: https://github.com/mcs07/MolVS/releases
.. _`MolVS Travis CI`: https://travis-ci.org/mcs07/MolVS
