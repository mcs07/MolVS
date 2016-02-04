#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for fragment.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem

from molvs.standardize import Standardizer
from molvs.fragment import FragmentRemover


logging.basicConfig(level=logging.DEBUG)


def fragment_parent_smiles(smiles, prefer_organic=False):
    """Utility function that returns the fragment parent SMILES for given a SMILES string."""
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Standardizer(prefer_organic=prefer_organic).fragment_parent(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def fragment_removal_smiles(smiles, leave_last=True):
    """Utility function that returns the result SMILES after FragmentRemover is applied to given a SMILES string."""
    mol = Chem.MolFromSmiles(smiles.encode('utf8'))
    mol = FragmentRemover(leave_last=leave_last).remove(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def test_fragment_parent():
    """Fragments created by standardization breaking metal-nonmetal covalent bonds."""
    assert fragment_parent_smiles('[Na]OC(=O)c1ccccc1') == 'O=C([O-])c1ccccc1'


def test_fragment_parent2():
    """Fragments created by standardization breaking metal-nonmetal covalent bonds."""
    assert fragment_parent_smiles('c1ccccc1C(=O)O[Ca]OC(=O)c1ccccc1') == 'O=C([O-])c1ccccc1'


def test_fragment_parent3():
    """Fragments created by standardization breaking metal-nonmetal covalent bonds."""
    assert fragment_parent_smiles('[Pt](Cl)(Cl)(O)(O)(NC(C)C)NC(C)C') == 'CC(C)[NH-]'


def test_fragment_parent4():
    """Mercury containing compound."""
    assert fragment_parent_smiles('CC[Hg]SC1=C(C=CC=C1)C(=O)[O][Na]') == 'CC[Hg]Sc1ccccc1C(=O)[O-]'


def test_fragment_parent5():
    """Covalent bond with metal."""
    assert fragment_parent_smiles('[Ag]OC(=O)O[Ag]') == 'O=C([O-])[O-]'


def test_fragment_parent6():
    """Salt without charges."""
    assert fragment_parent_smiles('[Na].O=C(O)c1ccccc1') == 'O=C([O-])c1ccccc1'


def test_fragment_parent7():
    """Multiple identical fragments."""
    assert fragment_parent_smiles('O=C(O)c1ccccc1.O=C(O)c1ccccc1.O=C(O)c1ccccc1') == 'O=C(O)c1ccccc1'


def test_fragment_parent8():
    """Multiple organic fragments of different sizes."""
    assert fragment_parent_smiles('O=C(O)CCC.O=C(O)CCCC.O=C(O)CCCCC.O=C(O)CCCC') == 'CCCCCC(=O)O'


def test_fragment_parent9():
    """No organic fragments."""
    assert fragment_parent_smiles('[N+](=O)([O-])[O-]') == 'O=[N+]([O-])[O-]'


def test_fragment_parent10():
    """No organic fragments."""
    assert fragment_parent_smiles('[N+](=O)([O-])[O-]', prefer_organic=True) == 'O=[N+]([O-])[O-]'


def test_fragment_parent11():
    """Larger inorganic fragment should be chosen."""
    assert fragment_parent_smiles('[N+](=O)([O-])[O-].[CH3+]') == 'O=[N+]([O-])[O-]'


def test_fragment_parent12():
    """Smaller organic fragment should be chosen over larger inorganic fragment."""
    assert fragment_parent_smiles('[N+](=O)([O-])[O-].[CH3+]', prefer_organic=True) == '[CH3+]'


def test_fragment_removal():
    """Single salt removal."""
    assert fragment_removal_smiles('CN(C)C.Cl') == 'CN(C)C'


def test_fragment_removal2():
    """Multiple salt removal."""
    assert fragment_removal_smiles('CN(C)C.Cl.Cl.Br') == 'CN(C)C'


def test_fragment_removal3():
    """FragmentPatterns should match entire fragments only, matches within larger fragments should be left."""
    assert fragment_removal_smiles('CN(Br)Cl') == 'CN(Cl)Br'


def test_fragment_removal4():
    """FragmentPatterns should match entire fragments only, matches within larger fragments should be left."""
    assert fragment_removal_smiles('CN(Br)Cl.Cl') == 'CN(Cl)Br'


def test_fragment_removal5():
    """Charged salts."""
    assert fragment_removal_smiles('C[NH+](C)(C).[Cl-]') == 'C[NH+](C)C'


def test_fragment_removal6():
    """Last match should be left."""
    assert fragment_removal_smiles('CC(=O)O.[Na]') == 'CC(=O)O'


def test_fragment_removal7():
    """Last match should be removed."""
    assert fragment_removal_smiles('CC(=O)O.[Na]', leave_last=False) == ''


def test_fragment_removal8():
    """Multiple identical last fragments should all be left."""
    assert fragment_removal_smiles('Cl.Cl') == 'Cl.Cl'


def test_fragment_removal9():
    """Test multiple fragment removal."""
    assert fragment_removal_smiles('[Na+].OC(=O)Cc1ccc(CN)cc1.OS(=O)(=O)C(F)(F)F') == 'NCc1ccc(CC(=O)O)cc1'


def test_fragment_removal10():
    """1,4-Dioxane should be removed."""
    assert fragment_removal_smiles('c1ccccc1O.O1CCOCC1') == 'Oc1ccccc1'


def test_fragment_removal11():
    """Benzene should be removed."""
    assert fragment_removal_smiles('c1ccccc1.CCCBr') == 'CCCBr'


def test_fragment_removal12():
    """Various fragments should be removed should be removed."""
    assert fragment_removal_smiles('CC(NC1=CC=C(O)C=C1)=O.CCCCC.O.CCO.CCCO.C1CCCCC1.C1CCCCCC1') == 'CC(=O)Nc1ccc(O)cc1'
