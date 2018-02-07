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


def test_reionize3():
    """"""
    mol = Chem.MolFromSmiles('C[N+]1=C[CH-]N(C(=N)N)/C1=C/[N+](=O)[O-]')
    r = Reionizer()
    mol = r.reionize(mol)
    assert Chem.MolToSmiles(mol) == 'C[N+]1=CCN(C(=N)N)C1=[C-][N+](=O)[O-]'


def test_should_complete():
    """Reionization should not infinitely loop forever on these molecules."""
    # GitHub Issue #14
    assert standardize_smiles('CCCCCCCCCCCCCCCCCC(=O)CC(=C)C(=O)O[Ti](=O)(OC(C)C)C(C)C') == 'C=C(CC(=O)[CH-]CCCCCCCCCCCCCCCC)C(=O)[O-].CC(C)[O-].CCC.[O-2].[Ti+5]'
    assert standardize_smiles('OP(=O)(O)[O-].OP(=O)([O-])[O-].[O-]S(=O)(=O)[O-].[Na+].[Na+].[Na+].[Mg+2].[Cl-].[Cl-].[K+].[K+]') == 'O=P([O-])(O)O.O=P([O-])([O-])O.O=S(=O)([O-])[O-].[Cl-].[Cl-].[K+].[K+].[Mg+2].[Na+].[Na+].[Na+]'


def test_forced_charge1():
    """Test forced charge correction maintaining overall neutral charge."""
    assert standardize_smiles('[Na].O=C(O)c1ccccc1') == 'O=C([O-])c1ccccc1.[Na+]'


def test_forced_charge2():
    """Test forced charge correction with no corresponding proton for neutralization."""
    # GitHub Issue #15
    assert standardize_smiles('[Na].[Na]') == '[Na+].[Na+]'
    # TODO: Arguably should become selenite ion... O=[Se]([O-])[O-]. Need an AcidBasePair?
    assert standardize_smiles('[Na].[Na].O[Se](O)=O') == 'O=[Se](O)O.[Na+].[Na+]'


# def test_reionize3():
#     """Test canonical ionization position when multiple equivalent possibilities."""
#     mol = Chem.MolFromSmiles('CC1=CC(=CC=C1S(O)=O)S([O-])=O')
#     mol2 = Chem.MolFromSmiles('CC1=CC(=CC=C1S([O-])=O)S(O)=O')
#     r = Reionizer()
#     mol = r.reionize(mol)
#     mol2 = r.reionize(mol2)
#     assert Chem.MolToSmiles(mol) == 'Cc1cc(S(=O)[O-])ccc1S(=O)O'
#     assert Chem.MolToSmiles(mol2) == 'Cc1cc(S(=O)[O-])ccc1S(=O)O'
#     assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(mol2)
#
#
# def test_reionize4():
#     """Test canonical ionization position when multiple equivalent possibilities."""
#     mol = Chem.MolFromSmiles('CCOC(=O)C(=O)[CH-]C#N')
#     mol2 = Chem.MolFromSmiles('[CH2-]COC(=O)C(=O)CC#N')
#     r = Reionizer()
#     mol = r.reionize(mol)
#     mol2 = r.reionize(mol2)
#     assert Chem.MolToSmiles(mol) == '[CH2-]COC(=O)C(=O)CC#N'
#     assert Chem.MolToSmiles(mol2) == ''
#     assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(mol2)
