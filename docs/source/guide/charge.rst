.. _charge:

Charges
=======

.. sectionauthor:: Matt Swain <m.swain@me.com>

This page gives details on dealing with charges in molecules.


Acid reionization
-----------------



Neutralization
--------------

- Get atoms:
    - `[+!H0!$(*~[-])]` positive with hydrogen
    - `[+H0!$(*~[-])]` quaternary positive
    - `[-!$(*~[+H0])]` negative total
    - `[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]` negative acid
- If quaternary positive charges
    - If more negative total than quaternary positive and negative acid > 0
        - Add hydrogen to first negative acid atom, increase formal charge
        - Repeat until quaternary positive ==negative total or no more negative acid
    -  TODO: Look into cases where acidic groups in zwitterions are inadvertently
       protonated
- else iterate over negative total atoms
    - Add hydrogens and increase formal change until neutral
- Iterate over positive with hydrogen atoms
    - Remove hydrogen and reduce formal change until neutral or no more hydrogens

