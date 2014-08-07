.. _api:

API documentation
=================

.. sectionauthor:: Matt Swain <m.swain@me.com>

.. module:: molvs

This part of the documentation is automatically generated from the MolVS source code and comments.

The MolVS package is made up of the following modules:

.. contents::
   :local:
   :depth: 1

.. automodule:: molvs.standardize
.. autoclass:: molvs.standardize.Standardizer(normalizations=NORMALIZATIONS, acid_base_pairs=ACID_BASE_PAIRS, tautomer_transforms=TAUTOMER_TRANSFORMS, tautomer_scores=TAUTOMER_SCORES, max_restarts=MAX_RESTARTS, max_tautomers=MAX_TAUTOMERS, prefer_organic=PREFER_ORGANIC)
   :special-members: __call__
   :members:
.. autofunction:: molvs.standardize.standardize_smiles
.. autofunction:: molvs.standardize.enumerate_tautomers_smiles
.. autofunction:: molvs.standardize.canonicalize_tautomer_smiles

.. automodule:: molvs.normalize
.. autodata:: molvs.normalize.NORMALIZATIONS
   :annotation:
.. autodata:: molvs.normalize.MAX_RESTARTS
   :annotation: = 200
.. autoclass:: molvs.normalize.Normalization
   :members:
.. autoclass:: molvs.normalize.Normalizer(normalizations=NORMALIZATIONS, max_restarts=MAX_RESTARTS)
   :special-members: __call__
   :members:

.. automodule:: molvs.metal
.. autoclass:: molvs.metal.MetalDisconnector
   :special-members: __call__
   :members:

.. automodule:: molvs.tautomer
.. autodata:: molvs.tautomer.TAUTOMER_TRANSFORMS
   :annotation:
.. autodata:: molvs.tautomer.TAUTOMER_SCORES
   :annotation:
.. autodata:: molvs.tautomer.MAX_TAUTOMERS
   :annotation: = 1000
.. autoclass:: molvs.tautomer.TautomerTransform
   :members:
.. autoclass:: molvs.tautomer.TautomerScore
   :members:
.. autoclass:: molvs.tautomer.TautomerCanonicalizer(transforms=TAUTOMER_TRANSFORMS, scores=TAUTOMER_SCORES, max_tautomers=MAX_TAUTOMERS)
   :special-members: __call__
   :members:
.. autoclass:: molvs.tautomer.TautomerEnumerator(transforms=TAUTOMER_TRANSFORMS, max_tautomers=MAX_TAUTOMERS)
   :special-members: __call__
   :members:

.. automodule:: molvs.fragment
.. autodata:: molvs.fragment.REMOVE_FRAGMENTS
   :annotation:
.. autodata:: molvs.fragment.LEAVE_LAST
   :annotation: = True
.. autodata:: molvs.fragment.PREFER_ORGANIC
   :annotation: = False
.. autoclass:: molvs.fragment.FragmentPattern
   :members:
.. autofunction:: molvs.fragment.is_organic
.. autoclass:: molvs.fragment.FragmentRemover(fragments=REMOVE_FRAGMENTS, leave_last=LEAVE_LAST)
   :special-members: __call__
   :members:
.. autoclass:: molvs.fragment.LargestFragmentChooser(prefer_organic=PREFER_ORGANIC)
   :special-members: __call__
   :members:

.. automodule:: molvs.charge
.. autodata:: molvs.charge.ACID_BASE_PAIRS
   :annotation:
.. autoclass:: molvs.charge.AcidBasePair
   :members:
.. autoclass:: molvs.charge.Reionizer(acid_base_pairs=ACID_BASE_PAIRS)
   :special-members: __call__
   :members:
.. autoclass:: molvs.charge.Uncharger
   :special-members: __call__
   :members:

.. automodule:: molvs.validate
.. autodata:: molvs.validate.SIMPLE_FORMAT
   :annotation: = '%(levelname)s: [%(validation)s] %(message)s'
.. autodata:: molvs.validate.LONG_FORMAT
   :annotation: = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'
.. autoclass:: molvs.validate.Validator(validations=VALIDATIONS, log_format=SIMPLE_FORMAT, level=logging.INFO, stdout=False, raw=False)
   :special-members: __call__
   :members:
.. autofunction:: molvs.validate.validate_smiles

.. _molvs_validations:
.. automodule:: molvs.validations
.. autodata:: molvs.validations.VALIDATIONS
   :annotation:
.. autoclass:: molvs.validations.Validation
   :members:
.. autoclass:: molvs.validations.SmartsValidation
   :members:
.. autoclass:: molvs.validations.IsNoneValidation
.. autoclass:: molvs.validations.NoAtomValidation
.. autoclass:: molvs.validations.DichloroethaneValidation
.. autoclass:: molvs.validations.FragmentValidation
.. autoclass:: molvs.validations.NeutralValidation
.. autoclass:: molvs.validations.IsotopeValidation

.. automodule:: molvs.cli

.. automodule:: molvs.errors
.. autoexception:: molvs.errors.MolVSError()
.. autoexception:: molvs.errors.StandardizeError()
.. autoexception:: molvs.errors.ValidateError()
.. autoexception:: molvs.errors.StopValidateError()
