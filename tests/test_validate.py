#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for charge.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

from molvs.validate import validate_smiles


def test_none():
    """IsNoneValidation should log due to SMILES parse error."""
    assert validate_smiles('3478q439g98h') == ['ERROR: [IsNoneValidation] Molecule is None']


def test_no_atoms():
    """An empty SMILES produces a mol with not atoms."""
    assert validate_smiles('') == [u'ERROR: [NoAtomValidation] No atoms are present']


def test_fragment():
    """FragmentValidation should identify 1,2-dichloroethane."""
    assert validate_smiles('ClCCCl.c1ccccc1O') == ['INFO: [FragmentValidation] 1,2-dichloroethane is present']


def test_fragment2():
    """FragmentValidation should identify 1,2-dimethoxyethane."""
    assert validate_smiles('COCCOC.CCCBr') == ['INFO: [FragmentValidation] 1,2-dimethoxyethane is present']


def test_charge():
    """NeutralValidation should identify net overall charge."""
    assert validate_smiles('O=C([O-])c1ccccc1') == ['INFO: [NeutralValidation] Not an overall neutral system (-1)']
    assert validate_smiles('CN=[NH+]CN=N') == ['INFO: [NeutralValidation] Not an overall neutral system (+1)']


def test_isotope():
    """IsotopeValidation should identify atoms with isotope labels."""
    assert validate_smiles('[13CH4]') == ['INFO: [IsotopeValidation] Molecule contains isotope 13C']
    assert validate_smiles('[2H]C(Cl)(Cl)Cl') == ['INFO: [IsotopeValidation] Molecule contains isotope 2H']
    assert validate_smiles('[2H]OC([2H])([2H])[2H]') == ['INFO: [IsotopeValidation] Molecule contains isotope 2H']
