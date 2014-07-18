.. MolVS documentation master file, created by sphinx-quickstart on Thu Apr 24 14:35:38 2014.

MolVS: Molecule Validation and Standardization
==============================================

.. sectionauthor:: Matt Swain <m.swain@me.com>

**MolVS** is a molecule validation and standardization tool, written in Python using the RDKit chemistry framework.

Building a collection of chemical structures from different sources can be difficult due to differing representations,
drawing conventions and mistakes. MolVS can standardize chemical structures to improve data quality, help with
de-duplication and identify relationships between molecules.

The available standardization tasks include disconnecting metals, normalizing functional groups, reionizing partially
ionized acids, discarding salts and solvents, choosing a canonical tautomer, neutralizing charges and more.

There are sensible defaults that make it easy to get started::

    >>> from molvs import standardize_smiles
    >>> standardize_smiles('C2(=C1C(=NC=N1)[NH]C(=N2)N)O')
    'Nc1nc(O)c2ncnc-2[nH]1'

Each standardization module is also available separately, allowing the development of custom standardization processes.

Features
--------

TODO

User guide
----------

A step-by-step guide to getting started with MolVS.

.. toctree::
   :maxdepth: 2

   guide/intro
   guide/install
   guide/gettingstarted
   guide/validate
   guide/standardize
   guide/tautomer
   guide/fragment
   guide/charge

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method. This is automatically generated
from the MolVS source code and comments.

.. toctree::
   :maxdepth: 2

   api


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`

