# -*- coding: utf-8 -*-
"""validate.py"""

from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import logging


# TODO: None of this works yet. Still figuring out exactly how to do this.


VALIDATIONS = (

)


class Validation(object):
    """The base class that all Validations must inherit from."""

    def apply(self, mol, logger):
        self.adapter = logging.LoggerAdapter(logger, {'validation': type(self).__name__})
        self.run(mol)

    def run(self, mol):
        """

        :return: The level for this Validation
        """
        raise NotImplementedError("Validation subclasses must implement the run method")



class TmpValidation(Validation):

    def run(self, mol):
        self.adapter.warn('A message from TmpValidation here')


# TODO: Think about logging formats
SIMPLE_FORMAT = '%(levelname)s: %(message)s'
LONG_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'


class LogHandler(logging.Handler):
    """A simple logging Handler that just stores logs in an array until flushed."""

    def __init__(self):
        logging.Handler.__init__(self)
        self.logs = []

    def emit(self, record):
        """Append the record."""
        self.logs.append(record)

    def flush(self):
        """Clear the log records."""
        self.acquire()
        try:
            self.buffer = []
        finally:
            self.release()

    def close(self):
        """Close the handler."""
        self.flush()
        logging.Handler.close(self)


class Validator(object):

    def __init__(self, validations=VALIDATIONS, log_format=SIMPLE_FORMAT):
        self.validations = validations
        self.log_format = log_format
        self.log = logging.getLogger(type(self).__name__)
        self.handler = LogHandler()
        self.handler.setFormatter(logging.Formatter(self.log_format))
        self.log.addHandler(self.handler)

    def validate(self, mol):
        self.handler.flush()
        self.log.warn('message here')
        self.log.info('message here2')
        self.log.warn('message here3')
        print(self.handler.logs)
        for log in self.handler.logs:
            print(log)




v = Validator()
v.validate(None)



# Attempt to read molecule using RDKit
# - Can produce errors, e.g. invalid atom valence
# - Is it possible to log the normalization of e.g. pentavalent nitro?

# - WARN/ERROR: Are all atoms defined/real - no query atoms or invalid elements, r-group things
# - WARN: No atoms present
# - INFO: Not an overall neutral system
# - INFO: Contains isotopes
# - INFO: Contains unknown stereo (Perform stereochemistry perception first?)
# - INFO: Nonstandard tautomer (log SMILES of tautomer parent, or the name of the tautomer transform?)
# - WARN: InChI generation failed
# - WARN: Contains covalent bond to metal (that would be broken by MetalDisconnector)
# - WARN: Contains solvent molecules (in addition other fragment)
# - WARN: More than 99 rings causes problems with SMILES
# - INFO: Cis azo dye is unusual
# - WARN: Adjacent atoms with like charges (i.e. both positive or both negative)
# - INFO: Has more than one radical centre
# - INFO: ethane, methane molecules present
# - INFO: Boron, Sulfur atoms with no explicit bonds
# - INFO: Solvent molecules present (only if also other fragments)
# - INFO: One unknown stereocentre and no defined stereocentres (probably racemate, so info not warn)
# - WARN: More than one undefined stereocentre and no defined stereocentres
# - INFO: One undefined stereocentre and at least one defined stereocentre (epimer or mixture of anomers, so info not warn)
# - WARN: More than one undefined stereocentre and at least one defined stereocentre
# - INFO: Unknown double bond stereochemistry
# - WARN: Ring containing stereobonds?
# - INFO: Not canonical tautomer


# Coordinates?
# Info - Lack of coordinates? Uneven bond lengths?

# Web services (needs to be optional)
# Info - Could not match to ChemSpider ID, PubChem CID
# UniChem from EBI could be useful here, otherwise use each API directly




# Allow definition of MolSchema to set custom validations on e.g.

# People can define a filterer
# This has a series of validations, and the required output - e.g. no error or no warns?

