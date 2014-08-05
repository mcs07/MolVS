#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for charge.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division

import nose
from nose.tools import eq_
from rdkit import Chem

from molvs.validate import validate_smiles, Validator


def test_none():
    """"""
    eq_(validate_smiles(''), ['ERROR: [IsNoneValidation] Molecule is None'])


def test_dichloroethane():
    """"""
    eq_(validate_smiles('ClCCCl.c1ccccc1O'), [u'INFO: [FragmentValidation] 1,2-dichloroethane is present'])


def test_dimethoxyethane():
    """"""
    eq_(validate_smiles('COCCOC.CCCBr'), [u'INFO: [FragmentValidation] 1,2-dimethoxyethane is present'])


if __name__ == '__main__':
    nose.main()
