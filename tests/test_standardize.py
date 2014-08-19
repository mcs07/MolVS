#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for standardize.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

import nose
from nose.tools import eq_

from molvs.standardize import standardize_smiles


logging.basicConfig(level=logging.DEBUG)


def test_aromaticity():
    """Check aromaticity is correctly perceived."""
    eq_(standardize_smiles('C1=CC=CC=C1'), 'c1ccccc1')


def test_aromaticity2():
    """Both rings should be aromatic."""
    eq_(standardize_smiles('C[N]1C=NC2=C1C(=O)N(C)C(=O)N2C'), 'Cn1cnc2c1c(=O)n(C)c(=O)n2C')
    eq_(standardize_smiles('Cn1cnc2c1c(=O)n(C)c(=O)n2C'), 'Cn1cnc2c1c(=O)n(C)c(=O)n2C')


def test_aromaticity3():
    """Redo incorrect aromatization."""
    eq_(standardize_smiles('C=Cc1ccc2c(c1)NC(=O)/C/2=C\\c1ccc[nH]1'), 'C=Cc1ccc2c(c1)NC(=O)/C2=C\\c1ccc[nH]1')
    eq_(standardize_smiles('C=Cc1ccc2c(c1)[nH]c(=O)/c/2=C\\c1ccc[nH]1'), 'C=Cc1ccc2c(c1)NC(=O)/C2=C\\c1ccc[nH]1')


def test_stereochemistry():
    """Check stereochemistry is correctly perceived."""
    eq_(standardize_smiles('Cl\\C=C/Cl'), 'Cl/C=C\\Cl')


def test_disconnect_metal():
    """Break metal-organic covalent bonds."""
    eq_(standardize_smiles('[Na]OC(=O)c1ccccc1'), '[Na+].O=C([O-])c1ccccc1')


def test_disconnect_metal2():
    """Break metal-organic covalent bonds."""
    eq_(standardize_smiles('c1ccccc1C(=O)O[Ca]OC(=O)c1ccccc1'), '[Ca+2].O=C([O-])c1ccccc1.O=C([O-])c1ccccc1')


def test_disconnect_metal3():
    """Disconnect Pt in metal complex."""
    eq_(standardize_smiles('[Pt](Cl)(Cl)(O)(O)(NC(C)C)NC(C)C'), '[OH-].[OH-].[Cl-].[Cl-].[Pt+6].CC(C)[NH-].CC(C)[NH-]')


def test_disconnect_metal4():
    """Leave mercury covalently bonded."""
    eq_(standardize_smiles('CC[Hg]SC1=C(C=CC=C1)C(=O)[O][Na]'), '[Na+].CC[Hg]Sc1ccccc1C(=O)[O-]')


def test_disconnect_metal5():
    """Silver carbonate. Unsure about this one."""
    eq_(standardize_smiles('[Ag]OC(=O)O[Ag]'), '[Ag+].[Ag+].O=C([O-])[O-]')


def test_charge_free_metal():
    """Charge free neutral metal with carboxylic acid."""
    eq_(standardize_smiles('[Na].O=C(O)c1ccccc1'), '[Na+].O=C([O-])c1ccccc1')


def test_nitro_normalization():
    """Normalize nitro group."""
    eq_(standardize_smiles('C1(=CC=CC=C1)[N+](=O)[O-]'), 'O=[N+]([O-])c1ccccc1')


def test_nitro_normalization2():
    """Normalize nitro group."""
    eq_(standardize_smiles('O=[N](=O)c1ccccc1'), 'O=[N+]([O-])c1ccccc1')


def test_nitro_normalization3():
    """Normalize nitro group."""
    eq_(standardize_smiles('[O-][N+](=O)c1ccccc1'), 'O=[N+]([O-])c1ccccc1')


def test_nitro_normalization4():
    """Normalize nitro group."""
    eq_(standardize_smiles('[N](=O)(=O)O'), 'O=[N+]([O-])O')


def test_nitro_normalization5():
    """Normalize nitro group."""
    eq_(standardize_smiles('O[N+](=O)[O-]'), 'O=[N+]([O-])O')


def test_pyridine_oxide_normalization():
    """Normalize pyridine oxide."""
    eq_(standardize_smiles('C1=[N](C=CC=C1)=O'), '[O-][n+]1ccccc1')


def test_pyridine_oxide_normalization2():
    """Normalize pyridine oxide."""
    eq_(standardize_smiles('O=n1ccccc1'), '[O-][n+]1ccccc1')


def test_sulfone_normalization():
    """Normalize sulfone."""
    eq_(standardize_smiles('C[S+2]([O-])([O-])C'), 'CS(C)(=O)=O')


def test_sulfone_normalization2():
    """Normalize sulfone."""
    eq_(standardize_smiles('C[S+2]([O-])([O-])O'), 'CS(=O)(=O)O')


def test_sulfoxide_normalization():
    """Normalize sulfoxide."""
    eq_(standardize_smiles('CS(=O)C'), 'C[S+](C)[O-]')


def test_sulfoxide_normalization2():
    """Normalize sulfoxide."""
    eq_(standardize_smiles('COC1=CC2=C(C=C1)[N]C(=N2)[S](=O)CC3=C(C(=C(C=N3)C)OC)C'), 'COc1ccc2[n]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1')


def test_sulfoxide_normalization3():
    """Normalize sulfoxide."""
    eq_(standardize_smiles('COc1ccc2c(c1)nc([nH]2)S(=O)Cc1ncc(c(c1C)OC)C'), 'COc1ccc2[nH]c([S+]([O-])Cc3ncc(C)c(OC)c3C)nc2c1')


def test_azide_normalization():
    """Normalize azide."""
    eq_(standardize_smiles('C1(=CC=C(C=C1)N)N=[N]#N'), '[N-]=[N+]=Nc1ccc(N)cc1')


def test_diazo_normalization():
    """Normalize diazo."""
    eq_(standardize_smiles('[N](#N)=C1C(NC(N=C1)=O)=O'), '[N-]=[N+]=C1C=NC(=O)NC1=O')


def test_phosphate_normalization():
    """Normalize phosphate."""
    eq_(standardize_smiles('C1=NC=C([N]1)CO[P+]([O-])([O-])[O-]'), 'O=P([O-])([O-])OCc1cnc[n]1')


def test_hydrazine_diazonium_normalization():
    """Normalize hydrazine-diazonium."""
    eq_(standardize_smiles('CNNC[N+]#N'), 'CN=[NH+]CN=N')


def test_amidinium_normalization():
    """Normalize amidinium."""
    eq_(standardize_smiles('[C+](C)(N)N'), 'CC(N)=[NH2+]')


def test_multi_fragment_normalization():
    """All fragments should stay if one gets transformed by normalization."""
    eq_(standardize_smiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1'), '[Na+].O=C([O-])c1ccc(CS(=O)=O)cc1')
    eq_(standardize_smiles('[Na+].[O-]C(=O)c1ccc(C[S+2]([O-])([O-]))cc1'), '[Na+].O=C([O-])c1ccc(CS(=O)=O)cc1')


def test_1_3_nonaromatic_charge_recombination():
    """Recombine non-aromatic 1,3-separated charges."""
    eq_(standardize_smiles('C[N-]C(C)=[N+](C)C'), 'CN=C(C)N(C)C')


def test_1_3_aromatic_charge_recombination():
    """Recombine aromatic 1,3-separated charges."""
    eq_(standardize_smiles('[n-]1c(=[N+](C)C)cccc1'), 'CN(C)c1ccccn1')


def test_1_3_aromatic_charge_recombination2():
    """Recombine aromatic 1,3-separated charges."""
    eq_(standardize_smiles('C[n+]1c([N-](C))cccc1'), 'CN=c1ccccn1C')


def test_pyrimidone_charge_recombination():
    """Recombine aromatic 1,3-separated charges to form pyrimidone."""
    eq_(standardize_smiles('[O-]c1[n+](C)cccc1'), 'Cn1ccccc1=O')


def test_pyrimidone_charge_recombination2():
    """Recombine aromatic 1,3-separated charges to form pyrimidone."""
    eq_(standardize_smiles('COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c([O-])[n+](C)c3c2cc1OC'), 'COc1cc2ccc3c4c(OC)cc(OC)c(OC)c4c(=O)n(C)c3c2cc1OC')


def test_1_5_nonaromatic_charge_recombination():
    """Recombine non-aromatic 1,5-separated charges."""
    eq_(standardize_smiles('C[N-]C=CC=[N+](C)C'), 'CN=CC=CN(C)C')


def test_1_5_aromatic_charge_recombination():
    """Recombine aromatic 1,5-separated charges."""
    eq_(standardize_smiles('[n-]1ccc(=[N+](C)C)cc1'), 'CN(C)c1ccncc1')


def test_1_5_aromatic_charge_recombination2():
    """Recombine aromatic 1,5-separated charges."""
    eq_(standardize_smiles('C[n+]1ccc([N-]C)cc1'), 'CN=c1ccn(C)cc1')


def test_charge_to_protonated_atom():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('CNC=[N+](C)C'), 'C[NH+]=CN(C)C')


def test_charge_to_protonated_atom2():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('[nH]1c(=[N+](C)C)cccc1'), 'CN(C)c1cccc[nH+]1')


def test_charge_to_protonated_atom3():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('CNc1[n+](C)cccc1'), 'C[NH+]=c1ccccn1C')


def test_charge_to_protonated_atom4():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('CNc1[n+](C)cco1'), 'C[NH+]=c1occn1C')


def test_charge_to_protonated_atom5():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('CNC=CC=[N+](C)C'), 'C[NH+]=CC=CN(C)C')


def test_charge_to_protonated_atom6():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('[nH]1ccc(=[N+](C)C)cc1'), 'CN(C)c1cc[nH+]cc1')


def test_charge_to_protonated_atom7():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('CNc1cc[n+](C)cc1'), 'C[NH+]=c1ccn(C)cc1')


def test_charge_to_protonated_atom8():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('C[n+]1ccc2[nH]ccc2c1'), 'Cn1ccc2[nH+]ccc-2c1')


def test_charge_to_protonated_atom9():
    """Shift positive charge from nonprotonated to protonated atom."""
    eq_(standardize_smiles('CNC=CC=[N+](C)C'), 'C[NH+]=CC=CN(C)C')


def test_transform_maintains_ring():
    """Ensure no transforms inadvertently breaks open rings."""
    eq_(standardize_smiles('[nH]1ccc2cccc[n+]12'), 'c1cc2ccccn2[nH+]1')


def test_equal_reionize():
    """Don't change partially ionized acid with two equally strong acid groups."""
    eq_(standardize_smiles('C1=C(C=CC(=C1)[S]([O-])(=O)=O)[S](O)(=O)=O'), 'O=S(=O)([O-])c1ccc(S(=O)(=O)O)cc1')


def test_reionize():
    """Partially ionized acid where proton should be moved to weaker acid."""
    eq_(standardize_smiles('C1=C(C=CC(=C1)[S]([O-])=O)[S](O)(=O)=O'), 'O=S(O)c1ccc(S(=O)(=O)[O-])cc1')


def test_reionize2():
    """Partially ionized acid where proton should be moved to weaker acid."""
    eq_(standardize_smiles('C1=C(C=CC(=C1)[P]([O-])(=O)O)[S](O)(=O)=O'), 'O=P(O)(O)c1ccc(S(=O)(=O)[O-])cc1')


def test_reionize3():
    """Partially ionized acid where proton should be moved to weaker acid."""
    eq_(standardize_smiles('C1=C(C=C(C(=C1)O)C(=O)[O-])[S](O)(=O)=O.[Na+]'), '[Na+].O=C(O)c1cc(S(=O)(=O)[O-])ccc1O')


def test_charge_preservation():
    """Unusually charged histidine should be preserved. (See charge parent for normalization)."""
    eq_(standardize_smiles('[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]'), '[NH3+]C(Cc1cnc[nH]1)C(=O)[O-]')


def test_charge_preservation2():
    """"""
    eq_(standardize_smiles('[Cl-].C[NH+](C)C'), '[Cl-].C[NH+](C)C')


def test_disconnect_metal6():
    """"""
    eq_(standardize_smiles('C1(CCCCC1)[Zn]Br'), '[Zn+2].[Br-].[CH-]1CCCCC1')


if __name__ == '__main__':
    nose.main()
