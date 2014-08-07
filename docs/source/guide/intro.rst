.. _introduction:

Introduction
============

.. sectionauthor:: Matt Swain <m.swain@me.com>

Building a collection of chemical structures from various different sources is difficult. There are differing
file formats, molecular representations, drawing conventions, and things that are just plain wrong.

A lot of this arises due to our chemical models being an imperfect description of reality, but even within the idealized
models there is often no single correct answer to whether two differently represented molecules are actually "the same".
Whether tautomers or isomers of the same molecule should be considered equivalent or distinct entities can depend
entirely on the specific application.

MolVS tries to address this problem through customizable validation and standardization processes, combined with the
concept of "parent" molecule relationships to allow multiple simultaneous degrees of standardization.

This guide provides a quick tour through MolVS concepts and functionality.


MolVS license
-------------

MolVS is released under the MIT License. This is a short, permissive software license that allows commercial use,
modifications, distribution, sublicensing and private use. Basically, you can do whatever you want with MolVS as long as
you include the original copyright and license in any copies or derivative projects.

See the `LICENSE file`_ for the full text of the license.

.. _`LICENSE file`: https://github.com/mcs07/MolVS/blob/master/LICENSE
