#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for metal.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem

from molvs.metal import MetalDisconnector
from molvs.standardize import standardize_smiles


logging.basicConfig(level=logging.DEBUG)


def test_standardize():
    """Test covalent metal is disconnected during standardize."""
    assert standardize_smiles('CCC(=O)O[Na]') == 'CCC(=O)[O-].[Na+]'


def test_standardize2():
    """Test metal ion is untouched during standardize."""
    assert standardize_smiles('CCC(=O)[O-].[Na+]') == 'CCC(=O)[O-].[Na+]'


def test_standardize3():
    """Test Hg is disconnected from O during standardize."""
    assert standardize_smiles('CCC(=O)O[Hg]') == 'CCC(=O)[O-].[Hg+]'


def test_standardize4():
    """Test dimethylmercury is not disconnected during standardize."""
    assert standardize_smiles('C[Hg]C') == 'C[Hg]C'


def test_standardize5():
    """Test zirconium (IV) ethoxide."""
    assert standardize_smiles('CCO[Zr](OCC)(OCC)OCC') == 'CC[O-].CC[O-].CC[O-].CC[O-].[Zr+4]'


def test_standardize6():
    """Test Grignard reagent."""
    # TODO: Should we disconnect this?
    assert standardize_smiles('c1ccccc1[Mg]Br') == 'Br[Mg]c1ccccc1'


def test_metaldisconnector1():
    """Test direct usage of MetalDisconnector class."""
    mol = Chem.MolFromSmiles('NC(CC(=O)O)C(=O)[O-].O.O.[Na+]')
    md = MetalDisconnector()
    mol = md.disconnect(mol)
    assert Chem.MolToSmiles(mol) == 'NC(CC(=O)O)C(=O)[O-].O.O.[Na+]'


def test_metaldisconnector2():
    """Test direct usage of MetalDisconnector class."""
    mol = Chem.MolFromSmiles('CCC(=O)O[Na]')
    md = MetalDisconnector()
    mol = md.disconnect(mol)
    assert Chem.MolToSmiles(mol) == 'CCC(=O)[O-].[Na+]'
