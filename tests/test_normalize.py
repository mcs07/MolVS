#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for normalize.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem

from molvs.normalize import Normalizer


logging.basicConfig(level=logging.DEBUG)


def normalize_smiles(smiles):
    """Utility function that runs normalization rules on a given a SMILES string."""
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Normalizer().normalize(mol)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_nitro():
    """Test nitro group normalization. Actually handled automatically by RDKit."""
    assert normalize_smiles('CN(=O)=O') == 'C[N+](=O)[O-]'


def test_sulfoxide():
    """Test sulfoxide normalization."""
    assert normalize_smiles('CS(C)=O') == 'C[S+](C)[O-]'


def test_sulfone():
    """"""
    assert normalize_smiles('C[S+2]([O-])([O-])O') == 'CS(=O)(=O)O'


def test_1_3_charge_recombination():
    """Test 1,3-separated charges are recombined."""
    assert normalize_smiles('CC([O-])=[N+](C)C') == 'CC(=O)N(C)C'


def test_1_3_charge_recombination_aromatic():
    """Test 1,3-separated charges are recombined."""
    assert normalize_smiles('C[n+]1ccccc1[O-]') == 'Cn1ccccc1=O'


def test_1_3_charge_recombination_exception():
    """Test a case where 1,3-separated charges should not be recombined."""
    assert normalize_smiles('CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]') == 'CC12CCCCC1(Cl)[N+]([O-])=[N+]2[O-]'


def test_1_5_charge_recombination():
    """Test 1,5-separated charges are recombined."""
    assert normalize_smiles('C[N+](C)=C\C=C\[O-]') == 'CN(C)C=CC=O'


def test_1_5_charge_recombination_exception():
    """Test a case where 1,5-separated charges should not be recombined."""
    assert normalize_smiles('C[N+]1=C2C=[N+]([O-])C=CN2CCC1') == 'C[N+]1=C2C=[N+]([O-])C=CN2CCC1'
