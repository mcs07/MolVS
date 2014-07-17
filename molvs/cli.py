# -*- coding: utf-8 -*-
"""
molvs.cli
~~~~~~~~~

This module contains a command line interface for standardization.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import argparse
import logging

from molvs import standardize_smiles


log = logging.getLogger(__name__)


def main():
    """Main function for molvs command line interface."""
    # For now, this just accepts a SMILES input argument and returns a standardized SMILES to stdout
    # TODO: Add standardization options, file input/output (formats?), tautomer enumeration
    parser = argparse.ArgumentParser()
    parser.add_argument('-:', nargs='?', help='input smiles', dest='smiles')
    args = parser.parse_args()
    print(standardize_smiles(args.smiles))

# TODO: Consider separate standardize and validate command line tools
