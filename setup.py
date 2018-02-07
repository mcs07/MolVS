#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup


if os.path.exists('README.rst'):
    long_description = open('README.rst').read()
else:
    long_description = '''MolVS is a molecule validation and standardization tool, written in Python using the RDKit
chemistry framework. Building a collection of chemical structures from different sources can be difficult due to
differing representations, drawing conventions and mistakes. MolVS can standardize chemical structures to improve data
quality, help with de-duplication and identify relationships between molecules.
'''

setup(
    name='MolVS',
    version='0.1.0',
    author='Matt Swain',
    author_email='m.swain@me.com',
    license='MIT',
    url='https://github.com/mcs07/MolVS',
    packages=['molvs'],
    description='Molecule Validation and Standardization',
    long_description=long_description,
    keywords='chemistry cheminformatics rdkit',
    zip_safe=False,
    tests_require=['pytest'],
    install_requires=['six'],
    entry_points={'console_scripts': ['molvs = molvs.cli:main']},
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
