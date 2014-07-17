#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for charge.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

import nose
from nose.tools import eq_
from rdkit import Chem

from molvs.standardize import Standardizer


logging.basicConfig(level=logging.DEBUG)


def charge_parent_smiles(smiles, prefer_organic=False):
    """Utility function that returns the charge parent SMILES for given a SMILES string."""
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Standardizer(prefer_organic=prefer_organic).charge_parent(mol)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_charge_parent():
    """Test neutralization of ionized acids and bases."""
    eq_(charge_parent_smiles('C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])'), 'NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O')


def test_charge_parent2():
    """Test preservation of zwitterion."""
    eq_(charge_parent_smiles('n(C)1cc[n+]2cccc([O-])c12'), 'Cn1cc[n+]2cccc([O-])c12')


def test_charge_parent3():
    """Choline should be left with a positive charge."""
    eq_(charge_parent_smiles('C[N+](C)(C)CCO'), 'C[N+](C)(C)CCO')


def test_charge_parent4():
    """This should have the hydrogen removed to give deanol as a charge parent."""
    eq_(charge_parent_smiles('C[NH+](C)CCO'), 'CN(C)CCO')


def test_charge_parent5():
    """Sodium benzoate to benzoic acid."""
    eq_(charge_parent_smiles('[Na+].O=C([O-])c1ccccc1'), 'O=C(O)c1ccccc1')


def test_charge_parent6():
    """Benzoate ion to benzoic acid."""
    eq_(charge_parent_smiles('O=C([O-])c1ccccc1'), 'O=C(O)c1ccccc1')


def test_charge_parent7():
    """Charges in histidine should be neutralized."""
    eq_(charge_parent_smiles('[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]'), 'NC(Cc1cnc[nH]1)C(=O)O')


def test_charge_parent8():
    """"""
    eq_(charge_parent_smiles('C[NH+](C)(C).[Cl-]'), 'CN(C)C')


def test_charge_parent9():
    """No organic fragments."""
    eq_(charge_parent_smiles('[N+](=O)([O-])[O-]'), 'O=[N+]([O-])[O-]')



def test_charge_parent10():
    """No organic fragments."""
    eq_(charge_parent_smiles('[N+](=O)([O-])[O-]', prefer_organic=True), 'O=[N+]([O-])[O-]')


def test_charge_parent11():
    """Larger inorganic fragment should be chosen."""
    eq_(charge_parent_smiles('[N+](=O)([O-])[O-].[CH2]'), 'O=[N+]([O-])[O-]')


def test_charge_parent12():
    """Smaller organic fragment should be chosen over larger inorganic fragment."""
    eq_(charge_parent_smiles('[N+](=O)([O-])[O-].[CH2]', prefer_organic=True), '[CH2]')


if __name__ == '__main__':
    nose.main()
