.. _cli:

Command Line Tool
=================

.. sectionauthor:: Matt Swain <m.swain@me.com>

MolVS comes with a simple command line tool that allows standardization and validation by typing ``molvs`` at the
command line.

Standardization
---------------

See standardization help by typing ``molvs standardize -h``::

    usage: molvs standardize [infile] [-i {smi,mol,sdf}] [-O <outfile>]
                             [-o {smi,mol,sdf}] [-: <smiles>]

    positional arguments:
      infile                input filename

    optional arguments:
      -i {smi,mol,sdf}, --intype {smi,mol,sdf}
                            input filetype
      -: <smiles>, --smiles <smiles>
                            input SMILES instead of file
      -O <outfile>, --outfile <outfile>
                            output filename
      -o {smi,mol,sdf}, --outtype {smi,mol,sdf}
                            output filetype




Validation
----------

See validation help by typing ``molvs validate -h``::

    usage: molvs validate [infile] [-i {smi,mol,sdf}] [-O <outfile>]
                          [-: <smiles>]

    positional arguments:
      infile                input filename

    optional arguments:
      -i {smi,mol,sdf}, --intype {smi,mol,sdf}
                            input filetype
      -: <smiles>, --smiles <smiles>
                            input SMILES instead of file
      -O <outfile>, --outfile <outfile>
                            output filename


Examples
--------

SMILES standardization:

.. code-block:: bash

    $ molvs standardize -:"C[n+]1c([N-](C))cccc1"
    CN=c1ccccn1C

Specifying an output format:

.. code-block:: bash

    $ molvs standardize -:"[N](=O)(=O)O" -o mol

         RDKit

      4  3  0  0  0  0  0  0  0  0999 V2000
        0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      1  3  2  0
      1  4  1  0
    M  CHG  2   1   1   2  -1
    M  END

Using stdin:

.. code-block:: bash

    $ echo "C[n+]1c([N-](C))cccc1" | molvs standardize
    CN=c1ccccn1C

Specifying an input file:

.. code-block:: bash

    $ molvs standardize example.mol
    CN=c1ccccn1C

Specifying an output file:

.. code-block:: bash

    $ molvs standardize example.mol -O output.smi
    $ molvs standardize example.mol -O output.mol
    $ molvs standardize example.mol -O output -o mol

Logging validations to stdout:

.. code-block:: bash

    $ molvs validate -:"O=C([O-])c1ccccc1"
    INFO: [NeutralValidation] Not an overall neutral system (-1)

Logging validations to a file:

.. code-block:: bash

    $ molvs validate -:"O=C([O-])c1ccccc1" -O logs.txt
