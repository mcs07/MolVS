# -*- coding: utf-8 -*-
"""
molvs.fragment
~~~~~~~~~~~~~~

This module contains tools for dealing with molecules with more than one covalently bonded unit. The main classes are
:class:`~molvs.fragment.LargestFragmentChooser`, which returns the largest covalent unit in a molecule, and
:class:`~molvs.fragment.FragmentRemover`, which filters out fragments from a molecule using SMARTS patterns.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from molvs.utils import memoized_property


log = logging.getLogger(__name__)


class FragmentPattern(object):
    """A fragment defined by a SMARTS pattern."""

    def __init__(self, name, smarts):
        """Initialize a FragmentPattern with a name and a SMARTS pattern.

        :param name: A name for this FragmentPattern.
        :param smarts: A SMARTS pattern.
        """
        self.name = name
        self.smarts_str = smarts

    @memoized_property
    def smarts(self):
        return Chem.MolFromSmarts(self.smarts_str.encode('utf8'))

    def __repr__(self):
        return 'FragmentPattern({!r}, {!r})'.format(self.name, self.smarts_str)

    def __str__(self):
        return self.name


#: The default list of :class:`FragmentPatterns <molvs.fragment.FragmentPattern>` to be used by
#: :class:`~molvs.fragment.FragmentRemover`.
REMOVE_FRAGMENTS = (
    FragmentPattern('Fluorine', '[F]'),
    FragmentPattern('Chlorine', '[Cl]'),
    FragmentPattern('Bromine', '[Br]'),
    FragmentPattern('Iodine', '[I]'),
    FragmentPattern('Lithium', '[Li]'),
    FragmentPattern('Sodium', '[Na]'),
    FragmentPattern('Potassium', '[K]'),
    FragmentPattern('Calcium', '[Ca]'),
    FragmentPattern('Magnesium', '[Mg]'),
    FragmentPattern('Aluminium', '[Al]'),
    FragmentPattern('Barium', '[Ba]'),
    FragmentPattern('Bismuth', '[Bi]'),
    FragmentPattern('Silver', '[Ag]'),
    FragmentPattern('Strontium', '[Sr]'),
    FragmentPattern('Zinc', '[Zn]'),
    FragmentPattern('Ammonia/Ammonium', 'N'),
    FragmentPattern('Water/Hydroxide', 'O'),
    FragmentPattern('Methyl amine', 'CN'),
    FragmentPattern('Sulfide', 'S'),
    FragmentPattern('Nitrate', '[N](=O)(O)O'),
    FragmentPattern('Phosphate', '[P](=O)(O)(O)O'),
    FragmentPattern('Hexafluorophosphate', '[P](F)(F)(F)(F)(F)F'),
    FragmentPattern('Sulfate', '[S](=O)(=O)(O)O'),
    FragmentPattern('Methyl sulfonate', '[CH3][S](=O)(=O)(O)'),
    FragmentPattern('p-toluene sulfonate', 'c1cc([CH3])ccc1[S](=O)(=O)(O)'),
    FragmentPattern('trifluorosulfonate', 'OS(=O)(=O)C(F)(F)F'),
    FragmentPattern('Acetate', '[CH3]C(=O)O'),
    FragmentPattern('TFA', 'FC(F)(F)C(=O)O'),
    FragmentPattern('Fumarate/Maleate', 'OC(=O)C=CC(=O)O'),
    FragmentPattern('Oxalate', 'OC(=O)C(=O)O'),
    FragmentPattern('Tartrate', 'OC(=O)C(O)C(O)C(=O)O'),
    FragmentPattern('Dicylcohexylammonium', 'C1CCCCC1[NH]C1CCCCC1'),
    FragmentPattern('Ethanol', 'CCO'),
    FragmentPattern('Acetone', 'CC(=O)C'),
    FragmentPattern('DMSO', 'CS(=O)C'),
    #FragmentPattern('', ''),
    # TODO: This list is far from complete
)

#: The default value for whether to ensure at least one fragment is left after FragmentRemover is applied.
LEAVE_LAST = True

#: The default value for whether LargestFragmentChooser sees organic fragments as "larger" than inorganic fragments.
PREFER_ORGANIC = False


def is_organic(fragment):
    """Return true if fragment contains at least one carbon atom.

    :param fragment: The fragment as an RDKit Mol object.
    """
    # TODO: Consider a different definition?
    # Could allow only H, C, N, O, S, P, F, Cl, Br, I
    for a in fragment.GetAtoms():
        if a.GetAtomicNum() == 6:
            return True
    return False


class FragmentRemover(object):
    """A class for filtering out fragments using SMARTS patterns."""

    def __init__(self, fragments=REMOVE_FRAGMENTS, leave_last=LEAVE_LAST):
        """Initialize a FragmentRemover with an optional custom list of :class:`~molvs.fragment.FragmentPattern`.

        Setting leave_last to True will ensure at least one fragment is left in the molecule, even if it is matched by a
        :class:`~molvs.fragment.FragmentPattern`. Fragments are removed in the order specified in the list, so place
        those you would prefer to be left towards the end of the list. If all the remaining fragments match the same
        :class:`~molvs.fragment.FragmentPattern`, they will all be left.

        :param fragments: A list of :class:`~molvs.fragment.FragmentPattern` to remove.
        :param bool leave_last: Whether to ensure at least one fragment is left.
        """
        log.debug('Initializing FragmentRemover')
        self.fragments = fragments
        self.leave_last = leave_last

    def __call__(self, mol):
        """Calling a FragmentRemover instance like a function is the same as calling its remove(mol) method."""
        return self.remove(mol)

    def remove(self, mol):
        """Return the molecule with specified fragments removed.

        :param mol: The molecule to remove fragments from.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The molecule with fragments removed.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running FragmentRemover')
        # Iterate FragmentPatterns and remove matching fragments
        for frag in self.fragments:
            # If nothing is left or leave_last and only one fragment, end here
            if mol.GetNumAtoms() == 0 or (self.leave_last and len(Chem.GetMolFrags(mol)) <= 1):
                break
            # Apply removal for this FragmentPattern
            removed = Chem.DeleteSubstructs(mol, frag.smarts, onlyFrags=True)
            if self.leave_last and removed.GetNumAtoms() == 0:
                # All the remaining fragments match this pattern - leave them all
                break
            mol = removed
        return mol


class LargestFragmentChooser(object):
    """A class for selecting the largest covalent unit in a molecule with multiple fragments."""

    def __init__(self, prefer_organic=PREFER_ORGANIC):
        """

        If prefer_organic is set to True, any organic fragment will be considered larger than any inorganic fragment. A
        fragment is considered organic if it contains a carbon atom.

        :param bool prefer_organic: Whether to prioritize organic fragments above all others.
        """
        log.debug('Initializing LargestFragmentChooser')
        self.prefer_organic = prefer_organic

    def __call__(self, mol):
        """Calling a LargestFragmentChooser instance like a function is the same as calling its choose(mol) method."""
        return self.choose(mol)

    def choose(self, mol):
        """Return the largest covalent unit.

        The largest fragment is determined by number of atoms (including hydrogens). Ties are broken by taking the
        fragment with the higher molecular weight, and then by taking the first alphabetically by SMILES if needed.

        :param mol: The molecule to choose the largest fragment from.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The largest fragment.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running LargestFragmentChooser')
        # TODO: Alternatively allow a list of fragments to be passed as the mol parameter
        fragments = Chem.GetMolFrags(mol, asMols=True)
        largest = None
        for f in fragments:
            smiles = Chem.MolToSmiles(f, isomericSmiles=True)
            log.debug('Fragment: %s', smiles)
            organic = is_organic(f)
            if self.prefer_organic:
                # Skip this fragment if not organic and we already have an organic fragment as the largest so far
                if largest and largest['organic'] and not organic:
                    continue
                # Reset largest if it wasn't organic and this fragment is organic
                if largest and organic and not largest['organic']:
                    largest = None
            # Count atoms
            atoms = 0
            for a in f.GetAtoms():
                atoms += 1 + a.GetTotalNumHs()
            # Skip this fragment if fewer atoms than the largest
            if largest and atoms < largest['atoms']:
                continue
            # Skip this fragment if equal number of atoms but weight is lower
            weight = rdMolDescriptors.CalcExactMolWt(f)
            if largest and atoms == largest['atoms'] and weight < largest['weight']:
                continue
            # Skip this fragment if equal atoms and equal weight but smiles comes last alphabetically
            if largest and atoms == largest['atoms'] and weight == largest['weight'] and smiles > largest['smiles']:
                continue
            # Otherwise this is the largest so far
            log.debug('New largest fragment: %s (%s)', smiles, atoms)
            largest = {'smiles': smiles, 'fragment': f, 'atoms': atoms, 'weight': weight, 'organic': organic}
        return largest['fragment']
