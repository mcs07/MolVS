# -*- coding: utf-8 -*-
"""validate.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging
import sys

from rdkit import Chem

from .errors import StopValidateError
from .validations import VALIDATIONS


#: The default format for log messages.
SIMPLE_FORMAT = '%(levelname)s: [%(validation)s] %(message)s'

#: A more detailed format for log messages. Specify when initializing a Validator.
LONG_FORMAT = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'


class LogHandler(logging.Handler):
    """A simple logging Handler that just stores logs in an array until flushed."""

    def __init__(self):
        logging.Handler.__init__(self)
        self.logs = []

    @property
    def logmessages(self):
        return [self.format(record) for record in self.logs]

    def emit(self, record):
        """Append the record."""
        self.logs.append(record)

    def flush(self):
        """Clear the log records."""
        self.acquire()
        try:
            self.logs = []
        finally:
            self.release()

    def close(self):
        """Close the handler."""
        self.flush()
        logging.Handler.close(self)


class Validator(object):
    """"""

    def __init__(self, validations=VALIDATIONS, log_format=SIMPLE_FORMAT, level=logging.INFO, stdout=False, raw=False):
        """"""
        self.raw = raw
        # Set up logger and add default LogHandler
        self.log = logging.getLogger(type(self).__name__)
        self.log.setLevel(level)
        self.handler = LogHandler()
        self.handler.setFormatter(logging.Formatter(log_format))
        self.log.addHandler(self.handler)
        # Add stdout StreamHandler if specified in parameters
        if stdout:
            strhdlr = logging.StreamHandler(sys.stdout)
            strhdlr.setFormatter(logging.Formatter(log_format))
            self.log.addHandler(strhdlr)
        # Instantiate the validations
        self.validations = [validation(self.log) for validation in validations]

    def validate(self, mol):
        """"""
        # Clear any log messages from previous runs
        self.handler.flush()
        # Run every validation, stopping if StopValidateError is raised
        for validation in self.validations:
            try:
                validation(mol)
            except StopValidateError:
                break
        return self.handler.logs if self.raw else self.handler.logmessages


def validate_smiles(smiles):
    """Return log messages for a given SMILES string using the default validations.

    Note: This is a convenience function for quickly validating a single SMILES string. It is more efficient to use
    the :class:`~molvs.validate.Validator` class directly when working with many molecules or when custom options
    are needed.

    :param string smiles: The SMILES for the molecule.
    :returns: The SMILES for the standardized molecule.
    :rtype: list of strings.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles.encode('utf8'))
    logs = Validator().validate(mol)
    return logs
