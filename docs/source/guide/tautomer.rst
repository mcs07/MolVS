.. _tautomer:

Tautomers
=========

.. sectionauthor:: Matt Swain <m.swain@me.com>

This page gives details on tautomer enumeration and canonicalization.


Tautomer enumeration
--------------------



Tautomer canonicalization
-------------------------

- Enumerate all possible tautomers using transform rules.
- Use scoring system to determine canonical tautomer.
- After tautomer canonicalization, need to clean up stereochemistry:
- Remove stereochemistry from double bonds that are single in at least 1 tautomer
- TODO: Remove tetrahedral stereochemistry from sp3 atoms that are sp2 in at least 1 tautomer? (thalidomide?)
