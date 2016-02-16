#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for resonance.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem

from molvs.resonance import ResonanceEnumerator, enumerate_resonance_smiles


logging.basicConfig(level=logging.DEBUG)


def test_thiocyanate_ion():
    """"""
    assert enumerate_resonance_smiles('[S-]C#N') == {'N#C[S-]', '[N-]=C=S'}


def test_thiocyanate_ion2():
    """"""
    assert enumerate_resonance_smiles('[N-]=C=S') == {'N#C[S-]', '[N-]=C=S'}


def test_carbamimidoylbenzoic_acid():
    """Custom ResonanceEnumerate options allow unconstrained charges."""
    mol = Chem.MolFromSmiles('NC(=[NH2+])c1ccc(cc1)C(=O)[O-]')
    rs = ResonanceEnumerator().enumerate(mol)
    assert {Chem.MolToSmiles(r, isomericSmiles=True) for r in rs} == {'NC(=[NH2+])c1ccc(C(=O)[O-])cc1'}
    rs = ResonanceEnumerator(allow_incomplete_octets=True, unconstrained_anions=True, unconstrained_cations=True).enumerate(mol)
    assert len(rs) == 32


def test_mobile_charge():
    """"""
    assert enumerate_resonance_smiles('CN1CC[N+]2=C1C1=C(C=CC=C1)C1=CC=CC=C21') == {'CN1CC[n+]2c1c1ccccc1c1ccccc12', 'C[N+]1=C2c3ccccc3-c3ccccc3N2CC1'}


def test_mobile_charge2():
    """"""
    assert enumerate_resonance_smiles('C[N+]1=C2N(CC1)C1=CC=CC=C1C1=C2C=CC=C1') == {'CN1CC[n+]2c1c1ccccc1c1ccccc12', 'C[N+]1=C2c3ccccc3-c3ccccc3N2CC1'}




# Normalization limits: From p36 of the InChI technical manual
# If passed to standardizer, the local symmetry means that either N could get +ve charge, regardless of O and S locations?
# CN(C)C1=[NH+]C2=CC(=S)C(=O)C=C2N1
# C[N+](C)=C1NC2=CC(=O)C(=S)C=C2N1
# Can this be solved by resonance/tautomer enumeration and canonicalization?



