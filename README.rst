MolVS: Molecule Validation and Standardization
==============================================

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

Installation
------------

To install MolVS, simply run::

    pip install molvs

Alternatively, try one of the other `installation options`_.

Documentation
-------------

Full documentation is available at http://molvs.readthedocs.org.

Contribute
----------

-  Feature ideas and bug reports are welcome on the `Issue Tracker`_.
-  Fork the `source code`_ on GitHub, make changes and send a pull request.

License
-------

MolVS is licensed under the `MIT license`_.

Similar projects
----------------

There are a number of other projects with similar goals:

- Francis Atkinson's `Standardiser`_.
- The RSC's `Chemistry Validation and Standardization Platform (CVSP)`_.
- The `PubChem Standardization Service`_.
- Tripod `Structure standardizer`_.
- The `FDA Substance Registration System Standard Operating Procedure`_.


.. _`installation options`: http://molvs.readthedocs.org/en/latest/guide/install.html
.. _`source code`: https://github.com/mcs07/MolVS
.. _`Issue Tracker`: https://github.com/mcs07/MolVS/issues
.. _`MIT license`: https://github.com/mcs07/MolVS/blob/master/LICENSE
.. _`Standardiser`: https://wwwdev.ebi.ac.uk/chembl/extra/francis/standardiser/
.. _`Chemistry Validation and Standardization Platform (CVSP)`: http://cvsp.chemspider.com
.. _`PubChem Standardization Service`: https://pubchem.ncbi.nlm.nih.gov/standardize/standardize.cgi
.. _`Structure standardizer`: https://tripod.nih.gov/?p=61
.. _`FDA Substance Registration System Standard Operating Procedure`: http://www.fda.gov/downloads/ForIndustry/DataStandards/SubstanceRegistrationSystem-UniqueIngredientIdentifierUNII/ucm127743.pdf
