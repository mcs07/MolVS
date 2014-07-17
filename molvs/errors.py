# -*- coding: utf-8 -*-
"""
molvs.errors
~~~~~~~~~~~~

This module contains exceptions that are raised by MolVS.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""


class MolVSError(Exception):
    pass


class StandardizeError(MolVSError):
    pass


class ValidateError(MolVSError):
    pass
