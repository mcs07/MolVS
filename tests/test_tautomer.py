#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for tautomer.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

import nose
from nose.tools import eq_
from rdkit import Chem

from molvs.standardize import Standardizer, enumerate_tautomers_smiles, canonicalize_tautomer_smiles
from molvs.tautomer import TautomerEnumerator, TautomerCanonicalizer


logging.basicConfig(level=logging.DEBUG)


def test_1_3_keto_enol_enumeration():
    """Enumerate 1,3 keto/enol tautomer."""
    eq_(enumerate_tautomers_smiles('C1(=CCCCC1)O'), {'OC1=CCCCC1', 'O=C1CCCCC1'})


def test_1_3_keto_enol_enumeration2():
    """Enumerate 1,3 keto/enol tautomer."""
    eq_(enumerate_tautomers_smiles('C1(CCCCC1)=O'), {'OC1=CCCCC1', 'O=C1CCCCC1'})


def test_acetophenone_keto_enol_enumeration():
    """Enumerate acetophenone keto/enol tautomer."""
    eq_(enumerate_tautomers_smiles('C(=C)(O)C1=CC=CC=C1'), {'C=C(O)c1ccccc1', 'CC(=O)c1ccccc1'})


def test_acetone_keto_enol_enumeration2():
    """Enumerate acetone keto/enol tautomer."""
    eq_(enumerate_tautomers_smiles('CC(C)=O'), {'CC(C)=O', 'C=C(C)O'})


def test_keto_enol_enumeration():
    """keto/enol tautomer"""
    eq_(enumerate_tautomers_smiles('OC(C)=C(C)C'), {'C=C(O)C(C)C', 'CC(C)=C(C)O', 'CC(=O)C(C)C'})


def test_phenylpropanone_keto_enol_enumeration():
    """1-phenyl-2-propanone enol/keto"""
    eq_(enumerate_tautomers_smiles('c1(ccccc1)CC(=O)C'), {'C=C(O)Cc1ccccc1', 'CC(=O)Cc1ccccc1', 'CC(O)=Cc1ccccc1'})


def test_1_5_keto_enol_enumeration():
    """1,5 keto/enol tautomer"""
    eq_(enumerate_tautomers_smiles('Oc1nccc2cc[nH]c(=N)c12'), {'Nc1nccc2ccnc(O)c12', 'N=c1nccc2cc[nH]c(O)c1-2', 'Nc1[nH]ccc2ccnc(=O)c1-2', 'N=c1[nH]ccc2ccnc(O)c21', 'Nc1nccc2cc[nH]c(=O)c12', 'N=c1[nH]ccc2cc[nH]c(=O)c21'})


def test_1_5_keto_enol_enumeration2():
    """1,5 keto/enol tautomer"""
    eq_(enumerate_tautomers_smiles('C1(C=CCCC1)=O'), {'O=C1C=CCCC1', 'OC1=CCC=CC1', 'OC1=CC=CCC1', 'O=C1CC=CCC1', 'OC1=CCCC=C1'})


def test_1_5_keto_enol_enumeration3():
    """1,5 keto/enol tautomer"""
    eq_(enumerate_tautomers_smiles('C1(=CC=CCC1)O'), {'O=C1C=CCCC1', 'OC1=CCC=CC1', 'OC1=CC=CCC1', 'O=C1CC=CCC1', 'OC1=CCCC=C1'})


def test_aliphatic_imine_enumeration():
    """aliphatic imine tautomer"""
    eq_(enumerate_tautomers_smiles('C1(CCCCC1)=N'), {'N=C1CCCCC1', 'NC1=CCCCC1'})


def test_aliphatic_imine_enumeration2():
    """aliphatic imine tautomer"""
    eq_(enumerate_tautomers_smiles('C1(=CCCCC1)N'), {'N=C1CCCCC1', 'NC1=CCCCC1'})


def test_special_imine_enumeration():
    """special imine tautomer"""
    eq_(enumerate_tautomers_smiles('C1(C=CC=CN1)=CC'), {'CC=C1C=CC=CN1', 'CCc1ccccn1', 'CC=C1C=CCC=N1'})


def test_special_imine_enumeration2():
    """special imine tautomer"""
    eq_(enumerate_tautomers_smiles('C1(=NC=CC=C1)CC'), {'CC=C1C=CC=CN1', 'CCc1ccccn1', 'CC=C1C=CCC=N1'})


def test_1_3_aromatic_heteroatom_enumeration():
    """1,3 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('O=c1cccc[nH]1'), {'Oc1ccccn1', 'O=c1cccc[nH]1'})


def test_1_3_aromatic_heteroatom_enumeration2():
    """1,3 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1ccccn1'), {'Oc1ccccn1', 'O=c1cccc[nH]1'})


def test_1_3_aromatic_heteroatom_enumeration3():
    """1,3 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1ncc[nH]1'), {'Oc1ncc[nH]1', 'O=c1[nH]cc[nH]1'})


def test_1_3_heteroatom_enumeration():
    """1,3 heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('OC(C)=NC'), {'CN=C(C)O', 'CNC(C)=O', 'C=C(O)NC'})


def test_1_3_heteroatom_enumeration2():
    """1,3 heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('CNC(C)=O'), {'CN=C(C)O', 'CNC(C)=O', 'C=C(O)NC'})


def test_1_3_heteroatom_enumeration3():
    """1,3 heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('S=C(N)N'), {'N=C(N)S', 'NC(N)=S'})


def test_1_3_heteroatom_enumeration4():
    """1,3 heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('SC(N)=N'), {'N=C(N)S', 'NC(N)=S'})


def test_1_3_heteroatom_enumeration5():
    """1,3 heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('N=c1[nH]ccn(C)1'), {'Cn1ccnc1N', 'Cn1cc[nH]c1=N'})


def test_1_3_heteroatom_enumeration6():
    """1,3 heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('CN=c1[nH]cncc1'), {'CN=c1ccnc[nH]1', 'CNc1ccncn1', 'CN=c1cc[nH]cn1'})


def test_1_5_aromatic_heteroatom_enumeration():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1cccc2ccncc12'), {'O=c1cccc2cc[nH]cc1-2', 'Oc1cccc2ccncc12'})


def test_1_5_aromatic_heteroatom_enumeration2():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('O=c1cccc2cc[nH]cc1-2'), {'O=c1cccc2cc[nH]cc1-2', 'Oc1cccc2ccncc12'})


def test_1_5_aromatic_heteroatom_enumeration3():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Cc1n[nH]c2ncnn12'), {'C=C1NNc2ncnn21', 'Cc1n[nH]c2ncnn12', 'Cc1nnc2[nH]cnn12', 'C=C1NN=C2N=CNN12', 'Cc1nnc2nc[nH]n12', 'C=C1NN=C2NC=NN12'})


def test_1_5_aromatic_heteroatom_enumeration4():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Cc1nnc2nc[nH]n12'), {'C=C1NNc2ncnn21', 'Cc1n[nH]c2ncnn12', 'Cc1nnc2[nH]cnn12', 'C=C1NN=C2N=CNN12', 'Cc1nnc2nc[nH]n12', 'C=C1NN=C2NC=NN12'})


def test_1_5_aromatic_heteroatom_enumeration5():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1ccncc1'), {'Oc1ccncc1', 'O=c1cc[nH]cc1'})


def test_1_5_aromatic_heteroatom_enumeration6():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1c(cccc3)c3nc2ccncc12'), {'Oc1c2ccccc2nc2ccncc12', 'O=c1c2ccccc2nc2cc[nH]cc1-2', 'O=c1c2ccccc2[nH]c2ccncc21'})


def test_1_3_1_5_aromatic_heteroatom_enumeration():
    """1,3 and 1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1ncncc1'), {'Oc1ccncn1', 'O=c1ccnc[nH]1', 'O=c1cc[nH]cn1'})


def test_1_5_aromatic_heteroatom_enumeration7():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('C2(=C1C(=NC=N1)[NH]C(=N2)N)O'), {'N=c1[nH]c2ncnc-2c(O)[nH]1', 'Nc1nc2nc[nH]c2c(O)n1', 'N=c1nc(O)c2nc[nH]c2[nH]1', 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1', 'Nc1nc2ncnc-2c(O)[nH]1', 'N=c1nc2nc[nH]c2c(O)[nH]1', 'N=c1nc(O)c2[nH]cnc2[nH]1', 'Nc1nc(O)c2ncnc-2[nH]1', 'Nc1nc(=O)c2nc[nH]c2[nH]1', 'Nc1nc(=O)c2[nH]cnc2[nH]1', 'Nc1nc2[nH]cnc2c(O)n1', 'N=c1nc2[nH]cnc2c(O)[nH]1', 'Nc1nc2[nH]cnc2c(=O)[nH]1', 'Nc1nc2nc[nH]c2c(=O)[nH]1', 'N=c1[nH]c2nc[nH]c2c(=O)[nH]1'})


def test_1_5_aromatic_heteroatom_enumeration8():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O'), {'N=c1[nH]c2ncnc-2c(O)[nH]1', 'Nc1nc2nc[nH]c2c(O)n1', 'N=c1nc(O)c2nc[nH]c2[nH]1', 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1', 'Nc1nc2ncnc-2c(O)[nH]1', 'N=c1nc2nc[nH]c2c(O)[nH]1', 'N=c1nc(O)c2[nH]cnc2[nH]1', 'Nc1nc(O)c2ncnc-2[nH]1', 'Nc1nc(=O)c2nc[nH]c2[nH]1', 'Nc1nc(=O)c2[nH]cnc2[nH]1', 'Nc1nc2[nH]cnc2c(O)n1', 'N=c1nc2[nH]cnc2c(O)[nH]1', 'Nc1nc2[nH]cnc2c(=O)[nH]1', 'Nc1nc2nc[nH]c2c(=O)[nH]1', 'N=c1[nH]c2nc[nH]c2c(=O)[nH]1'})


def test_1_5_aromatic_heteroatom_enumeration9():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Oc1n(C)ncc1'), {'Cn1nccc1O', 'CN1N=CCC1=O', 'Cn1[nH]ccc1=O'})


def test_1_5_aromatic_heteroatom_enumeration10():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('O=c1nc2[nH]ccn2cc1'), {'O=c1ccn2cc[nH]c2n1', 'Oc1ccn2ccnc2n1', 'O=c1ccn2ccnc2[nH]1'})


def test_1_5_aromatic_heteroatom_enumeration11():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('N=c1nc[nH]cc1'), {'N=c1cc[nH]cn1', 'N=c1ccnc[nH]1', 'Nc1ccncn1'})


def test_1_5_aromatic_heteroatom_enumeration12():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('N=c(c1)ccn2cc[nH]c12'), {'N=c1ccn2cc[nH]c2c1', 'Nc1ccn2ccnc2c1'})


def test_1_5_aromatic_heteroatom_enumeration13():
    """1,5 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('CN=c1nc[nH]cc1'), {'CN=c1ccnc[nH]1', 'CNc1ccncn1', 'CN=c1cc[nH]cn1'})


def test_1_7_aromatic_heteroatom_enumeration():
    """1,7 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1'), {'c1ccc2c(c1)=NC(C1=NC3C=CC=CC3=N1)N=2', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1', 'c1ccc2[nH]c(C3=NC4C=CC=CC4=N3)nc2c1', 'c1ccc2[nH]c(C3N=c4ccccc4=N3)nc2c1', 'c1ccc2c(c1)=NC(=C1N=C3C=CC=CC3N1)N=2', 'c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2'})


def test_1_7_aromatic_heteroatom_enumeration2():
    """1,7 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2'), {'c1ccc2c(c1)=NC(C1=NC3C=CC=CC3=N1)N=2', 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1', 'c1ccc2[nH]c(C3=NC4C=CC=CC4=N3)nc2c1', 'c1ccc2[nH]c(C3N=c4ccccc4=N3)nc2c1', 'c1ccc2c(c1)=NC(=C1N=C3C=CC=CC3N1)N=2', 'c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2'})


def test_1_9_aromatic_heteroatom_enumeration():
    """1,9 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('CNc1ccnc2ncnn21'), {'CN=c1cc[nH]c2ncnn21', 'CN=c1ccnc2[nH]cnn21', 'CN=c1ccnc2nc[nH]n21', 'CNc1ccnc2ncnn21'})


def test_1_9_aromatic_heteroatom_enumeration2():
    """1,9 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('CN=c1ccnc2nc[nH]n21'), {'CN=c1cc[nH]c2ncnn21', 'CN=c1ccnc2[nH]cnn21', 'CN=c1ccnc2nc[nH]n21', 'CNc1ccnc2ncnn21'})


def test_1_11_aromatic_heteroatom_enumeration():
    """1,11 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('Nc1ccc(C=C2C=CC(=O)C=C2)cc1'), {'Nc1ccc(C=C2C=CC(=O)C=C2)cc1', 'N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1', 'N=C1C=CC(=Cc2ccc(O)cc2)C=C1', 'N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1'})


def test_1_11_aromatic_heteroatom_enumeration2():
    """1,11 aromatic heteroatom H shift"""
    eq_(enumerate_tautomers_smiles('N=C1C=CC(=Cc2ccc(O)cc2)C=C1'), {'Nc1ccc(C=C2C=CC(=O)C=C2)cc1', 'N=C1C=CC(=CC2C=CC(=O)C=C2)C=C1', 'N=C1C=CC(=Cc2ccc(O)cc2)C=C1', 'N=C1C=CC(C=C2C=CC(=O)C=C2)C=C1'})


def test_heterocyclic_enumeration():
    """heterocyclic tautomer"""
    eq_(enumerate_tautomers_smiles('n1ccc2ccc[nH]c12'), {'c1cc2cccnc2[nH]1', 'c1cc2ccc[nH]c-2n1'})


def test_heterocyclic_enumeration2():
    """heterocyclic tautomer"""
    eq_(enumerate_tautomers_smiles('c1cc(=O)[nH]c2nccn12'), {'O=c1ccn2cc[nH]c2n1', 'Oc1ccn2ccnc2n1', 'O=c1ccn2ccnc2[nH]1'})


def test_heterocyclic_enumeration3():
    """heterocyclic tautomer"""
    eq_(enumerate_tautomers_smiles('c1cnc2c[nH]ccc12'), {'c1cc2cc[nH]cc-2n1', 'c1cc2ccncc2[nH]1'})


def test_heterocyclic_enumeration4():
    """heterocyclic tautomer"""
    eq_(enumerate_tautomers_smiles('n1ccc2c[nH]ccc12'), {'c1cc2cnccc2[nH]1', 'c1cc2c[nH]ccc-2n1'})


def test_heterocyclic_enumeration5():
    """heterocyclic tautomer"""
    eq_(enumerate_tautomers_smiles('c1cnc2ccc[nH]c12'), {'c1cc2[nH]cccc-2n1', 'c1cc2ncccc2[nH]1'})


def test_furanone_enumeration():
    """furanone tautomer"""
    eq_(enumerate_tautomers_smiles('C1=CC=C(O1)O'), {'Oc1ccco1', 'O=C1CC=CO1'})


def test_furanone_enumeration2():
    """furanone tautomer"""
    eq_(enumerate_tautomers_smiles('O=C1CC=CO1'), {'Oc1ccco1', 'O=C1CC=CO1'})


def test_keten_ynol_enumeration():
    """keten/ynol tautomer"""
    eq_(enumerate_tautomers_smiles('CC=C=O'), {'CC=C=O', 'CC#CO'})


def test_keten_ynol_enumeration2():
    """keten/ynol tautomer"""
    eq_(enumerate_tautomers_smiles('CC#CO'), {'CC=C=O', 'CC#CO'})


def test_ionic_nitro_aci_nitro_enumeration():
    """ionic nitro/aci-nitro tautomer"""
    eq_(enumerate_tautomers_smiles('C([N+](=O)[O-])C'), {'CC[N+](=O)[O-]', 'CC=[N+]([O-])O'})


def test_ionic_nitro_aci_nitro_enumeration2():
    """ionic nitro/aci-nitro tautomer"""
    eq_(enumerate_tautomers_smiles('C(=[N+](O)[O-])C'), {'CC[N+](=O)[O-]', 'CC=[N+]([O-])O'})


def test_oxim_nitroso_enumeration():
    """oxim nitroso tautomer"""
    eq_(enumerate_tautomers_smiles('CC(C)=NO'), {'CC(C)N=O', 'CC(C)=NO', 'C=C(C)NO'})


def test_oxim_nitroso_enumeration2():
    """oxim nitroso tautomer"""
    eq_(enumerate_tautomers_smiles('CC(C)N=O'), {'CC(C)N=O', 'CC(C)=NO', 'C=C(C)NO'})


def test_oxim_nitroso_enumeration3():
    """oxim/nitroso tautomer via phenol"""
    eq_(enumerate_tautomers_smiles('O=Nc1ccc(O)cc1'), {'O=NC1C=CC(=O)C=C1', 'O=C1C=CC(=NO)C=C1', 'O=Nc1ccc(O)cc1'})


def test_oxim_nitroso_enumeration4():
    """oxim/nitroso tautomer via phenol"""
    eq_(enumerate_tautomers_smiles('O=C1C=CC(=NO)C=C1'), {'O=NC1C=CC(=O)C=C1', 'O=C1C=CC(=NO)C=C1', 'O=Nc1ccc(O)cc1'})


def test_cyano_iso_cyanic_acid_enumeration():
    """cyano/iso-cyanic acid tautomer"""
    eq_(enumerate_tautomers_smiles('C(#N)O'), {'N#CO', 'N=C=O'})


def test_cyano_iso_cyanic_acid_enumeration2():
    """cyano/iso-cyanic acid tautomer"""
    eq_(enumerate_tautomers_smiles('C(=N)=O'), {'N#CO', 'N=C=O'})


def test_formamidinesulfinic_acid_enumeration():
    """formamidinesulfinic acid tautomer"""
    eq_(enumerate_tautomers_smiles('[S](=O)(=O)C(N)N'), {'N=C(N)S(=O)O', 'NC(N)S(=O)=O'})


def test_formamidinesulfinic_acid_enumeration2():
    """formamidinesulfinic acid tautomer"""
    eq_(enumerate_tautomers_smiles('[S](=O)(O)C(=N)N'), {'N=C(N)S(=O)O', 'NC(N)S(=O)=O'})


def test_isocyanide_enumeration():
    """isocyanide tautomer"""
    eq_(enumerate_tautomers_smiles('C#N'), {'[C-]#[NH+]', 'C#N'})


def test_isocyanide_enumeration2():
    """isocyanide tautomer"""
    eq_(enumerate_tautomers_smiles('[C-]#[NH+]'), {'[C-]#[NH+]', 'C#N'})


def test_phosphonic_acid_enumeration():
    """phosphonic acid tautomer"""
    eq_(enumerate_tautomers_smiles('[PH](=O)(O)(O)'), {'OP(O)O', 'O=[PH](O)O'})


def test_phosphonic_acid_enumeration2():
    """phosphonic acid tautomer"""
    eq_(enumerate_tautomers_smiles('P(O)(O)O'), {'OP(O)O', 'O=[PH](O)O'})


def test_mobile_double_stereochemistry_enumeration():
    """Remove stereochemistry from mobile double bonds"""
    eq_(enumerate_tautomers_smiles('c1(ccccc1)/C=C(/O)\\C'), {'C=C(O)Cc1ccccc1', 'CC(O)=Cc1ccccc1', 'CC(=O)Cc1ccccc1'})


def test_mobile_double_stereochemistry_enumeration2():
    """Remove stereochemistry from mobile double bonds"""
    eq_(enumerate_tautomers_smiles('C/C=C/C(C)=O'), {'C=C(O)C=CC', 'C=CCC(=C)O', 'CC=CC(C)=O', 'C=CCC(C)=O', 'C=CC=C(C)O'})


def test_mobile_double_stereochemistry_enumeration3():
    """Remove stereochemistry from mobile double bonds"""
    eq_(enumerate_tautomers_smiles('C/C=C\C(C)=O'), {'C=C(O)C=CC', 'C=CCC(=C)O', 'CC=CC(C)=O', 'C=CCC(C)=O', 'C=CC=C(C)O'})


def test_1_3_keto_enol_canonicalization():
    """1,3 keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(=CCCCC1)O'), 'O=C1CCCCC1')


def test_1_3_keto_enol_canonicalization2():
    """1,3 keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(CCCCC1)=O'), 'O=C1CCCCC1')


def test_acetophenone_keto_enol_canonicalization():
    """Acetophenone keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('C(=C)(O)C1=CC=CC=C1'), 'CC(=O)c1ccccc1')


def test_acetone_keto_enol_canonicalization():
    """Acetone keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('CC(C)=O'), 'CC(C)=O')


def test_keto_enol_canonicalization():
    """keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('OC(C)=C(C)C'), 'CC(=O)C(C)C')


def test_phenylpropanone_keto_enol_canonicalization():
    """1-phenyl-2-propanone enol/keto"""
    eq_(canonicalize_tautomer_smiles('c1(ccccc1)CC(=O)C'), 'CC(=O)Cc1ccccc1')


def test_1_5_keto_enol_canonicalization():
    """1,5 keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('Oc1nccc2cc[nH]c(=N)c12'), 'N=c1[nH]ccc2cc[nH]c(=O)c21')


def test_1_5_keto_enol_canonicalization2():
    """1,5 keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(C=CCCC1)=O'), 'O=C1C=CCCC1')


def test_1_5_keto_enol_canonicalization3():
    """1,5 keto/enol tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(=CC=CCC1)O'), 'O=C1C=CCCC1')


def test_aliphatic_imine_canonicalization():
    """aliphatic imine tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(CCCCC1)=N'), 'N=C1CCCCC1')


def test_aliphatic_imine_canonicalization2():
    """aliphatic imine tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(=CCCCC1)N'), 'N=C1CCCCC1')


def test_special_imine_canonicalization():
    """special imine tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(C=CC=CN1)=CC'), 'CCc1ccccn1')


def test_special_imine_canonicalization2():
    """special imine tautomer"""
    eq_(canonicalize_tautomer_smiles('C1(=NC=CC=C1)CC'), 'CCc1ccccn1')


def test_1_3_aromatic_heteroatom_canonicalization():
    """1,3 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('O=c1cccc[nH]1'), 'O=c1cccc[nH]1')


def test_1_3_aromatic_heteroatom_canonicalization2():
    """1,3 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1ccccn1'), 'O=c1cccc[nH]1')


def test_1_3_aromatic_heteroatom_canonicalization3():
    """1,3 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1ncc[nH]1'), 'O=c1[nH]cc[nH]1')


def test_1_3_heteroatom_canonicalization():
    """1,3 heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('OC(C)=NC'), 'CNC(C)=O')


def test_1_3_heteroatom_canonicalization2():
    """1,3 heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('CNC(C)=O'), 'CNC(C)=O')


def test_1_3_heteroatom_canonicalization3():
    """1,3 heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('S=C(N)N'), 'NC(N)=S')


def test_1_3_heteroatom_canonicalization4():
    """1,3 heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('SC(N)=N'), 'NC(N)=S')


def test_1_3_heteroatom_canonicalization5():
    """1,3 heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('N=c1[nH]ccn(C)1'), 'Cn1cc[nH]c1=N')


def test_1_3_heteroatom_canonicalization6():
    """1,3 heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('CN=c1[nH]cncc1'), 'CN=c1cc[nH]cn1')


def test_1_5_aromatic_heteroatom_canonicalization():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1cccc2ccncc12'), 'Oc1cccc2ccncc12')


def test_1_5_aromatic_heteroatom_canonicalization2():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('O=c1cccc2cc[nH]cc1-2'), 'Oc1cccc2ccncc12')


def test_1_5_aromatic_heteroatom_canonicalization3():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Cc1n[nH]c2ncnn12'), 'Cc1n[nH]c2ncnn12')


def test_1_5_aromatic_heteroatom_canonicalization4():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Cc1nnc2nc[nH]n12'), 'Cc1n[nH]c2ncnn12')


def test_1_5_aromatic_heteroatom_canonicalization5():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1ccncc1'), 'O=c1cc[nH]cc1')


def test_1_5_aromatic_heteroatom_canonicalization6():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1c(cccc3)c3nc2ccncc12'), 'O=c1c2ccccc2[nH]c2ccncc21')


def test_1_3_1_5_aromatic_heteroatom_canonicalization():
    """1,3 and 1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1ncncc1'), 'O=c1cc[nH]cn1')


def test_1_5_aromatic_heteroatom_canonicalization7():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('C2(=C1C(=NC=N1)[NH]C(=N2)N)O'), 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1')


def test_1_5_aromatic_heteroatom_canonicalization8():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('C2(C1=C([NH]C=N1)[NH]C(=N2)N)=O'), 'N=c1[nH]c2[nH]cnc2c(=O)[nH]1')


def test_1_5_aromatic_heteroatom_canonicalization9():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Oc1n(C)ncc1'), 'Cn1[nH]ccc1=O')


def test_1_5_aromatic_heteroatom_canonicalization10():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('O=c1nc2[nH]ccn2cc1'), 'O=c1ccn2cc[nH]c2n1')


def test_1_5_aromatic_heteroatom_canonicalization11():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('N=c1nc[nH]cc1'), 'N=c1cc[nH]cn1')


def test_1_5_aromatic_heteroatom_canonicalization12():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('N=c(c1)ccn2cc[nH]c12'), 'N=c1ccn2cc[nH]c2c1')


def test_1_5_aromatic_heteroatom_canonicalization13():
    """1,5 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('CN=c1nc[nH]cc1'), 'CN=c1cc[nH]cn1')


def test_1_7_aromatic_heteroatom_canonicalization():
    """1,7 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1'), 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1')


def test_1_7_aromatic_heteroatom_canonicalization2():
    """1,7 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('c1ccc2c(c1)NC(=C1N=c3ccccc3=N1)N2'), 'c1ccc2[nH]c(-c3nc4ccccc4[nH]3)nc2c1')


def test_1_9_aromatic_heteroatom_canonicalization():
    """1,9 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('CNc1ccnc2ncnn21'), 'CN=c1cc[nH]c2ncnn21')


def test_1_9_aromatic_heteroatom_canonicalization2():
    """1,9 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('CN=c1ccnc2nc[nH]n21'), 'CN=c1cc[nH]c2ncnn21')


def test_1_11_aromatic_heteroatom_canonicalization():
    """1,11 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('Nc1ccc(C=C2C=CC(=O)C=C2)cc1'), 'Nc1ccc(C=C2C=CC(=O)C=C2)cc1')


def test_1_11_aromatic_heteroatom_canonicalization2():
    """1,11 aromatic heteroatom H shift"""
    eq_(canonicalize_tautomer_smiles('N=C1C=CC(=Cc2ccc(O)cc2)C=C1'), 'Nc1ccc(C=C2C=CC(=O)C=C2)cc1')


def test_heterocyclic_canonicalization():
    """heterocyclic tautomer"""
    eq_(canonicalize_tautomer_smiles('n1ccc2ccc[nH]c12'), 'c1cc2cccnc2[nH]1')


def test_heterocyclic_canonicalization2():
    """heterocyclic tautomer"""
    eq_(canonicalize_tautomer_smiles('c1cc(=O)[nH]c2nccn12'), 'O=c1ccn2cc[nH]c2n1')


def test_heterocyclic_canonicalization3():
    """heterocyclic tautomer"""
    eq_(canonicalize_tautomer_smiles('c1cnc2c[nH]ccc12'), 'c1cc2ccncc2[nH]1')


def test_heterocyclic_canonicalization4():
    """heterocyclic tautomer"""
    eq_(canonicalize_tautomer_smiles('n1ccc2c[nH]ccc12'), 'c1cc2cnccc2[nH]1')


def test_heterocyclic_canonicalization5():
    """heterocyclic tautomer"""
    eq_(canonicalize_tautomer_smiles('c1cnc2ccc[nH]c12'), 'c1cc2ncccc2[nH]1')


def test_furanone_canonicalization():
    """furanone tautomer"""
    eq_(canonicalize_tautomer_smiles('C1=CC=C(O1)O'), 'Oc1ccco1')


def test_furanone_canonicalization2():
    """furanone tautomer"""
    eq_(canonicalize_tautomer_smiles('O=C1CC=CO1'), 'Oc1ccco1')


def test_keten_ynol_canonicalization():
    """keten/ynol tautomer"""
    eq_(canonicalize_tautomer_smiles('CC=C=O'), 'CC=C=O')


def test_keten_ynol_canonicalization2():
    """keten/ynol tautomer"""
    eq_(canonicalize_tautomer_smiles('CC#CO'), 'CC=C=O')


def test_ionic_nitro_aci_nitro_canonicalization():
    """ionic nitro/aci-nitro tautomer"""
    eq_(canonicalize_tautomer_smiles('C([N+](=O)[O-])C'), 'CC[N+](=O)[O-]')


def test_ionic_nitro_aci_nitro_canonicalization2():
    """ionic nitro/aci-nitro tautomer"""
    eq_(canonicalize_tautomer_smiles('C(=[N+](O)[O-])C'), 'CC[N+](=O)[O-]')


def test_oxim_nitroso_canonicalization():
    """oxim nitroso tautomer"""
    eq_(canonicalize_tautomer_smiles('CC(C)=NO'), 'CC(C)=NO')


def test_oxim_nitroso_canonicalization2():
    """oxim nitroso tautomer"""
    eq_(canonicalize_tautomer_smiles('CC(C)N=O'), 'CC(C)=NO')


def test_oxim_nitroso_phenol_canonicalization():
    """oxim/nitroso tautomer via phenol"""
    eq_(canonicalize_tautomer_smiles('O=Nc1ccc(O)cc1'), 'O=Nc1ccc(O)cc1')


def test_oxim_nitroso_phenol_canonicalization2():
    """oxim/nitroso tautomer via phenol"""
    eq_(canonicalize_tautomer_smiles('O=C1C=CC(=NO)C=C1'), 'O=Nc1ccc(O)cc1')


def test_cyano_iso_cyanic_acid_canonicalization():
    """cyano/iso-cyanic acid tautomer"""
    eq_(canonicalize_tautomer_smiles('C(#N)O'), 'N=C=O')


def test_cyano_iso_cyanic_acid_canonicalization2():
    """cyano/iso-cyanic acid tautomer"""
    eq_(canonicalize_tautomer_smiles('C(=N)=O'), 'N=C=O')


def test_formamidinesulfinic_acid_canonicalization():
    """formamidinesulfinic acid tautomer"""
    eq_(canonicalize_tautomer_smiles('[S](=O)(=O)C(N)N'), 'N=C(N)S(=O)O')


def test_formamidinesulfinic_acid_canonicalization2():
    """formamidinesulfinic acid tautomer"""
    eq_(canonicalize_tautomer_smiles('[S](=O)(O)C(=N)N'), 'N=C(N)S(=O)O')


def test_isocyanide_canonicalization():
    """isocyanide tautomer"""
    eq_(canonicalize_tautomer_smiles('C#N'), 'C#N')


def test_isocyanide_canonicalization2():
    """isocyanide tautomer"""
    eq_(canonicalize_tautomer_smiles('[C-]#[NH+]'), 'C#N')


def test_phosphonic_acid_canonicalization():
    """phosphonic acid tautomer"""
    eq_(canonicalize_tautomer_smiles('[PH](=O)(O)(O)'), 'O=[PH](O)O')


def test_phosphonic_acid_canonicalization2():
    """phosphonic acid tautomer"""
    eq_(canonicalize_tautomer_smiles('P(O)(O)O'), 'O=[PH](O)O')


if __name__ == '__main__':
    nose.main()
