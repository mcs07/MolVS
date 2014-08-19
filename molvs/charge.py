# -*- coding: utf-8 -*-
"""
molvs.charge
~~~~~~~~~~~~

This module implements tools for manipulating charges on molecules. In particular, :class:`~molvs.charge.Reionizer`,
which competitively reionizes acids such that the strongest acids ionize first, and :class:`~molvs.charge.Uncharger`,
which attempts to neutralize ionized acids and bases on a molecule.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import copy
import logging

from rdkit import Chem

from .utils import memoized_property


log = logging.getLogger(__name__)


class AcidBasePair(object):
    """An acid and its conjugate base, defined by SMARTS.

    A strength-ordered list of AcidBasePairs can be used to ensure the strongest acids in a molecule ionize first.
    """

    def __init__(self, name, acid, base):
        """Initialize an AcidBasePair with the following parameters:

        :param string name: A name for this AcidBasePair.
        :param string acid: SMARTS pattern for the protonated acid.
        :param string base: SMARTS pattern for the conjugate ionized base.
        """
        log.debug('Initializing AcidBasePair: %s', name)
        self.name = name
        self.acid_str = acid
        self.base_str = base

    @memoized_property
    def acid(self):
        log.debug('Loading AcidBasePair acid: %s', self.name)
        return Chem.MolFromSmarts(self.acid_str.encode('utf8'))

    @memoized_property
    def base(self):
        log.debug('Loading AcidBasePair base: %s', self.name)
        return Chem.MolFromSmarts(self.base_str.encode('utf8'))

    def __repr__(self):
        return 'AcidBasePair({!r}, {!r}, {!r})'.format(self.name, self.acid_str, self.base_str)

    def __str__(self):
        return self.name


#: The default list of AcidBasePairs, sorted from strongest to weakest. This list is derived from the Food and Drug
#: Administration Substance Registration System Standard Operating Procedure guide.
ACID_BASE_PAIRS = (
    AcidBasePair('-OSO3H', 'OS(=O)(=O)[OH]', 'OS(=O)(=O)[O-]'),
    AcidBasePair('–SO3H', '[!O]S(=O)(=O)[OH]', '[!O]S(=O)(=O)[O-]'),
    AcidBasePair('-OSO2H', 'O[SD3](=O)[OH]', 'O[SD3](=O)[O-]'),
    AcidBasePair('-SO2H', '[!O][SD3](=O)[OH]', '[!O][SD3](=O)[O-]'),
    AcidBasePair('-OPO3H2', 'OP(=O)([OH])[OH]', 'OP(=O)([OH])[O-]'),
    AcidBasePair('-PO3H2', '[!O]P(=O)([OH])[OH]', '[!O]P(=O)([OH])[O-]'),
    AcidBasePair('-CO2H', 'C(=O)[OH]', 'C(=O)[O-]'),
    AcidBasePair('thiophenol', 'c[SH]', 'c[S-]'),
    AcidBasePair('(-OPO3H)-', 'OP(=O)([OH])[O-]', 'OP(=O)([O-])[O-]'),
    AcidBasePair('(-PO3H)-', '[!O]P(=O)([OH])[O-]', '[!O]P(=O)([O-])[O-]'),
    AcidBasePair('phthalimide', 'O=C2c1ccccc1C(=O)[NH]2', 'O=C2c1ccccc1C(=O)[N-]2'),
    AcidBasePair('CO3H (peracetyl)', 'C(=O)O[OH]', 'C(=O)O[O-]'),
    AcidBasePair('alpha-carbon-hydrogen-nitro group', 'O=N(O)[CH]', 'O=N(O)[C-]'),
    AcidBasePair('-SO2NH2', 'S(=O)(=O)[NH2]', 'S(=O)(=O)[NH-]'),
    AcidBasePair('-OBO2H2', 'OB([OH])[OH]', 'OB([OH])[O-]'),
    AcidBasePair('-BO2H2', '[!O]B([OH])[OH]', '[!O]B([OH])[O-]'),
    AcidBasePair('phenol', 'c[OH]', 'c[O-]'),
    AcidBasePair('SH (aliphatic)', 'C[SH]', 'C[S-]'),
    AcidBasePair('(-OBO2H)-', 'OB([OH])[O-]', 'OB([O-])[O-]'),
    AcidBasePair('(-BO2H)-', '[!O]B([OH])[O-]', '[!O]B([O-])[O-]'),
    AcidBasePair('cyclopentadiene', '[CH2]1C=CC=C1', '[C-]1C=CC=C1'),
    AcidBasePair('-CONH2', 'C(=O)[NH2]', 'C(=O)[NH-]'),
    AcidBasePair('imidazole', 'c1cnc[n]1', 'c1cnc[n-]1'),
    AcidBasePair('-OH', '[CX4][OH]', '[CX4][O-]'),
    AcidBasePair('alpha-carbon-hydrogen-keto group', 'O=C[CH]', 'O=C[C-]'),
    AcidBasePair('alpha-carbon-hydrogen-acetyl ester group', 'OC(=O)[CH]', 'OC(=O)[C-]'),
    AcidBasePair('sp carbon hydrogen', 'C#[CH]', 'C#[C-]'),
    AcidBasePair('alpha-carbon-hydrogen-sulfone group', 'CS(=O)(=O)C[CH]', 'CS(=O)(=O)C[C-]'),
    AcidBasePair('alpha-carbon-hydrogen-sulfoxide group', 'C[SD3](=O)C[CH]', 'C[SD3](=O)C[C-]'),
    AcidBasePair('-NH2', '[CX4][NH2]', '[CX4][NH-]'),
    AcidBasePair('benzyl hydrogen', 'c[CD4H]', 'c[CD3-]'),
    AcidBasePair('sp2-carbon hydrogen', '[CX3]=[CX3H]', '[CX3]=[CX2-]'),
    AcidBasePair('sp3-carbon hydrogen', '[CX4H]', '[CX3-]'),
)


class Reionizer(object):
    """A class to reionize a molecule such that the strongest acids ionize first."""

    def __init__(self, acid_base_pairs=ACID_BASE_PAIRS):
        """Initialize a Reionizer with the following parameter:

        :param acid_base_pairs: A list of :class:`AcidBasePairs <molvs.charge.AcidBasePair>` to reionize, sorted from
                                strongest to weakest.
        """
        log.debug('Initializing Reionizer')
        self.acid_base_pairs = acid_base_pairs

    def __call__(self, mol):
        """Calling a Reionizer instance like a function is the same as calling its reionize(mol) method."""
        return self.reionize(mol)

    def reionize(self, mol):
        """If molecule with multiple acid groups is partially ionized, ensure strongest acids ionize first.

        The algorithm works as follows:

        - Use SMARTS to find the strongest protonated acid and the weakest ionized acid.
        - If the ionized acid is weaker than the protonated acid, swap proton and repeat.

        :param mol: The molecule to reionize.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The reionized molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running Reionizer')
        while True:
            ppos, poccur = self._strongest_protonated(mol)
            ipos, ioccur = self._weakest_ionized(mol)
            if ioccur and poccur and ppos < ipos:
                log.info('Moved proton from %s to %s', self.acid_base_pairs[ppos].name, self.acid_base_pairs[ipos].name)
                patom = mol.GetAtomWithIdx(poccur[-1])
                patom.SetFormalCharge(patom.GetFormalCharge() - 1)
                patom.SetNumExplicitHs(max(0, patom.GetNumExplicitHs() - 1))
                iatom = mol.GetAtomWithIdx(ioccur[-1])
                iatom.SetFormalCharge(iatom.GetFormalCharge() + 1)
                iatom.SetNumExplicitHs(iatom.GetTotalNumHs() + 1)
            else:
                Chem.SanitizeMol(mol)
                return mol

    def _strongest_protonated(self, mol):
        for position, pair in enumerate(self.acid_base_pairs):
            for occurrence in mol.GetSubstructMatches(pair.acid):
                return position, occurrence
        return None, None

    def _weakest_ionized(self, mol):
        for position, pair in enumerate(reversed(self.acid_base_pairs)):
            for occurrence in mol.GetSubstructMatches(pair.base):
                return len(self.acid_base_pairs) - position - 1, occurrence
        return None, None


class Uncharger(object):
    """Class for neutralizing ionized acids and bases.

    This class uncharges molecules by adding and/or removing hydrogens. For zwitterions, hydrogens are moved to
    eliminate charges where possible. However, in cases where there is a positive charge that is not neutralizable, an
    attempt is made to also preserve the corresponding negative charge.

    The method is derived from the neutralise module in `Francis Atkinson's standardiser tool
    <https://github.com/flatkinson/standardiser>`_, which is released under the Apache License v2.0.
    """

    def __init__(self):
        log.debug('Initializing Uncharger')
        #: Neutralizable positive charge (with hydrogens attached)
        self._pos_h = Chem.MolFromSmarts('[+!H0!$(*~[-])]'.encode('utf8'))
        #: Non-neutralizable positive charge (no hydrogens attached)
        self._pos_quat = Chem.MolFromSmarts('[+H0!$(*~[-])]'.encode('utf8'))
        #: Negative charge, not bonded to a positive charge with no hydrogens
        self._neg = Chem.MolFromSmarts('[-!$(*~[+H0])]'.encode('utf8'))
        #: Negative oxygen bonded to [C,P,S]=O, negative aromatic nitrogen?
        self._neg_acid = Chem.MolFromSmarts('[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]'.encode('utf8'))

    def __call__(self, mol):
        """Calling an Uncharger instance like a function is the same as calling its uncharge(mol) method."""
        return self.uncharge(mol)

    def uncharge(self, mol):
        """Neutralize molecule by adding/removing hydrogens. Attempts to preserve zwitterions.

        :param mol: The molecule to uncharge.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The uncharged molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running Uncharger')
        mol = copy.deepcopy(mol)
        # Get atom ids for matches
        p = [x[0] for x in mol.GetSubstructMatches(self._pos_h)]
        q = [x[0] for x in mol.GetSubstructMatches(self._pos_quat)]
        n = [x[0] for x in mol.GetSubstructMatches(self._neg)]
        a = [x[0] for x in mol.GetSubstructMatches(self._neg_acid)]
        # Neutralize negative charges
        if q:
            # Surplus negative charges more than non-neutralizable positive charges
            neg_surplus = len(n) - len(q)
            if a and neg_surplus > 0:
                # zwitterion with more negative charges than quaternary positive centres
                while neg_surplus > 0 and a:
                    # Add hydrogen to first negative acid atom, increase formal charge
                    # Until quaternary positive == negative total or no more negative acid
                    atom = mol.GetAtomWithIdx(a.pop(0))
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                    neg_surplus -= 1
                    log.info('Removed negative charge')
        else:
            #
            for atom in [mol.GetAtomWithIdx(x) for x in n]:
                while atom.GetFormalCharge() < 0:
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                    log.info('Removed negative charge')
        # Neutralize positive charges
        for atom in [mol.GetAtomWithIdx(x) for x in p]:
            # Remove hydrogen and reduce formal change until neutral or no more hydrogens
            while atom.GetFormalCharge() > 0 and atom.GetNumExplicitHs() > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                log.info('Removed positive charge')
        return mol
