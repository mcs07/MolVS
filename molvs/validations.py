# -*- coding: utf-8 -*-
"""
molvs.validations
~~~~~~~~~~~~~~~~~



:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem

from .errors import StopValidateError
from .fragment import FragmentPattern


class Validation(object):
    """The base class that all Validations must inherit from."""

    def __init__(self, log):
        self.log = logging.LoggerAdapter(log, {'validation': type(self).__name__})

    def __call__(self, mol):
        try:
            self.run(mol)
        except Exception as e:
            if isinstance(e, StopValidateError):
                raise e
            else:
                self.log.error('Validation failed: %s', e)

    def run(self, mol):
        """"""
        raise NotImplementedError("Validation subclasses must implement the run method")


class SmartsValidation(Validation):
    """Abstract class for Validations that log a message if a SMARTS pattern matches the molecule."""

    #: The logging level of the message.
    level = logging.INFO

    #: The message to log if the SMARTS pattern matches the molecule.
    message = 'Molecule matched %(smarts)s'

    #: Whether the SMARTS pattern should match an entire covalent unit.
    entire_fragment = False

    def __init__(self, log):
        super(SmartsValidation, self).__init__(log)
        self._smarts = Chem.MolFromSmarts(self.smarts)

    @property
    def smarts(self):
        """The SMARTS pattern as a string."""
        raise NotImplementedError('SmartsValidation subclasses must have a smarts attribute')

    def _check_matches(self, mol):
        if mol.HasSubstructMatch(self._smarts):
            self.log.log(self.level, self.message, {'smarts': self.smarts})

    def _check_matches_fragment(self, mol):
        matches = frozenset(frozenset(match) for match in mol.GetSubstructMatches(self._smarts))
        fragments = frozenset(frozenset(frag) for frag in Chem.GetMolFrags(mol))
        if matches & fragments:
            self.log.log(self.level, self.message, {'smarts': self.smarts})

    def run(self, mol):
        self.log.debug('Running %s', type(self).__name__)
        if self.entire_fragment:
            self._check_matches_fragment(mol)
        else:
            self._check_matches(mol)


class IsNoneValidation(Validation):
    """Validation to check whether None is passed to the Validator.

    This can happen if RDKit failed to parse an input format. If the molecule is None, no subsequent validations will
    run.
    """

    def run(self, mol):
        self.log.debug('Running IsNoneValidation')
        if mol is None:
            self.log.error('Molecule is None')
            raise StopValidateError()


class NoAtomValidation(Validation):
    """Validation to check whether the molecule has zero atoms.

    If the molecule has no atoms, no subsequent validations will run.
    """

    def run(self, mol):
        self.log.debug('Running NoAtomValidation')
        if mol.GetNumAtoms() == 0:
            self.log.warning('No atoms are present')
            raise StopValidateError()


class DichloroethaneValidation(SmartsValidation):
    """Validation to check if 1,2-Dichloroethane is present."""
    level = logging.INFO
    smarts = '[Cl]-[#6]-[#6]-[Cl]'
    entire_fragment = True
    message = '1,2-Dichloroethane is present'


class DimethoxyethaneValidation(SmartsValidation):
    """Validation to check if 1,2-Dimethoxyethane is present."""
    level = logging.INFO
    smarts = '[#6]-[#8]-[#6]-[#6]-[#8]-[#6]'
    entire_fragment = True
    message = '1,2-Dimethoxyethane is present'


VALIDATIONS = (
    IsNoneValidation,
    NoAtomValidation,
    DichloroethaneValidation,
    DimethoxyethaneValidation,
)



# - WARN/ERROR: Are all atoms defined/real - no query atoms or invalid elements, r-group things
# - INFO: Not an overall neutral system
# - INFO: Contains isotopes
# - INFO: Contains unknown stereo (Perform stereochemistry perception first?)
# - INFO: Nonstandard tautomer (log SMILES of tautomer parent, or the name of the tautomer transform?)
# - WARN: InChI generation failed
# - WARN: Contains covalent bond to metal (that would be broken by MetalDisconnector)
# - WARN: Contains solvent molecules (in addition other fragment)
# - WARN: More than 99 rings causes problems with SMILES
# - INFO: Cis azo dye is unusual
# - WARN: Adjacent atoms with like charges (i.e. both positive or both negative)
# - INFO: Has more than one radical centre
# - INFO: ethane, methane molecules present
# - INFO: Boron, Sulfur atoms with no explicit bonds
# - INFO: Solvent molecules present (only if also other fragments)
# - INFO: One unknown stereocentre and no defined stereocentres (probably racemate, so info not warn)
# - WARN: More than one undefined stereocentre and no defined stereocentres
# - INFO: One undefined stereocentre and at least one defined stereocentre (epimer or mixture of anomers, so info not warn)
# - WARN: More than one undefined stereocentre and at least one defined stereocentre
# - INFO: Unknown double bond stereochemistry
# - WARN: Ring containing stereobonds?
# - INFO: Not canonical tautomer


# Coordinates?
# Info - Lack of coordinates? Uneven bond lengths?

# Web services (needs to be optional)
# Info - Could not match to ChemSpider ID, PubChem CID
# UniChem from EBI could be useful here, otherwise use each API directly




# Allow definition of MolSchema to set custom validations on e.g.

# People can define a filterer
# This has a series of validations, and the required output - e.g. no error or no warns?

