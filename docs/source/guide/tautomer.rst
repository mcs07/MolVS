.. _tautomer:

Tautomers
=========

.. sectionauthor:: Matt Swain <m.swain@me.com>

This page gives details on tautomer enumeration and canonicalization.

Background
----------

Tautomers are sets of molecules that readily interconvert with each other through the movement of a hydrogen atom.
Tautomers have the same molecular formula and net charge, but they differ in terms of the positions of hydrogens and
the associated changes in adjacent double and single bonds.

Because they rapidly interconvert, for many applications tautomers are considered to be the same chemical compound.
And even in situations where it is important to treat tautomers as distinct compounds, it is still useful to be
aware of the tautomerism relationships between molecules in a collection.

Varying tautomeric forms of the same molecule can have significantly different fingerprints and descriptors, which can
negatively impact models for things like property prediction if they are used inconsistently.

There are two main tautomerism tasks that MolVS carries out:

- Tautomer enumeration: Finding the set of all the different possible tautomeric forms of a molecule.
- Tautomer canonicalization: Consistently picking one of the tautomers to be the canonical tautomer for the set.

Tautomer enumeration
--------------------

- All possible tautomers are generated using a series of transform rules.
- Remove stereochemistry from double bonds that are single in at least 1 tautomer.

Tautomer canonicalization
-------------------------

- Enumerate all possible tautomers using transform rules.
- Use scoring system to determine canonical tautomer.
- Canonical tautomer should be "reasonable" from a chemist's point of view, but isn't guaranteed to be the most energetically favourable.

