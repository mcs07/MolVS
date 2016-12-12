#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for charge.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem

from molvs.standardize import Standardizer, standardize_smiles
from molvs.charge import Reionizer


logging.basicConfig(level=logging.DEBUG)


def charge_parent_smiles(smiles, prefer_organic=False):
    """Utility function that returns the charge parent SMILES for given a SMILES string."""
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Standardizer(prefer_organic=prefer_organic).charge_parent(mol)
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_charge_parent():
    """Test neutralization of ionized acids and bases."""
    assert charge_parent_smiles('C(C(=O)[O-])(Cc1n[n-]nn1)(C[NH3+])(C[N+](=O)[O-])') == 'NCC(Cc1nn[nH]n1)(C[N+](=O)[O-])C(=O)O'


def test_charge_parent2():
    """Test preservation of zwitterion."""
    assert charge_parent_smiles('n(C)1cc[n+]2cccc([O-])c12') == 'Cn1cc[n+]2cccc([O-])c12'


def test_charge_parent3():
    """Choline should be left with a positive charge."""
    assert charge_parent_smiles('C[N+](C)(C)CCO') == 'C[N+](C)(C)CCO'


def test_charge_parent4():
    """This should have the hydrogen removed to give deanol as a charge parent."""
    assert charge_parent_smiles('C[NH+](C)CCO') == 'CN(C)CCO'


def test_charge_parent5():
    """Sodium benzoate to benzoic acid."""
    assert charge_parent_smiles('[Na+].O=C([O-])c1ccccc1') == 'O=C(O)c1ccccc1'


def test_charge_parent6():
    """Benzoate ion to benzoic acid."""
    assert charge_parent_smiles('O=C([O-])c1ccccc1') == 'O=C(O)c1ccccc1'


def test_charge_parent7():
    """Charges in histidine should be neutralized."""
    assert charge_parent_smiles('[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]') == 'NC(Cc1cnc[nH]1)C(=O)O'


def test_charge_parent8():
    """"""
    assert charge_parent_smiles('C[NH+](C)(C).[Cl-]') == 'CN(C)C'


def test_charge_parent9():
    """No organic fragments."""
    assert charge_parent_smiles('[N+](=O)([O-])[O-]') == 'O=[N+]([O-])[O-]'


def test_charge_parent10():
    """No organic fragments."""
    assert charge_parent_smiles('[N+](=O)([O-])[O-]', prefer_organic=True) == 'O=[N+]([O-])[O-]'


def test_charge_parent11():
    """Larger inorganic fragment should be chosen."""
    assert charge_parent_smiles('[N+](=O)([O-])[O-].[CH2]') == 'O=[N+]([O-])[O-]'


def test_charge_parent12():
    """Smaller organic fragment should be chosen over larger inorganic fragment."""
    assert charge_parent_smiles('[N+](=O)([O-])[O-].[CH2]', prefer_organic=True) == '[CH2]'


def test_standardize():
    """Test table salt."""
    assert standardize_smiles('[Na].[Cl]') == '[Cl-].[Na+]'


def test_reionize():
    """Test reionizer moves proton to weaker acid."""
    mol = Chem.MolFromSmiles('C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O')
    r = Reionizer()
    mol = r.reionize(mol)
    assert Chem.MolToSmiles(mol) == 'O=S(O)c1ccc(S(=O)(=O)[O-])cc1'


def test_reionize2():
    """Test charged carbon doesn't get recognised as alpha-carbon-hydrogen-keto."""
    mol = Chem.MolFromSmiles('CCOC(=O)C(=O)[CH-]C#N')
    r = Reionizer()
    mol = r.reionize(mol)
    assert Chem.MolToSmiles(mol) == 'CCOC(=O)C(=O)[CH-]C#N'
