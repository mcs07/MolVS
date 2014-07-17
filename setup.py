#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from setuptools import setup

import molvs


if os.path.exists('README.rst'):
    long_description = open('README.rst').read()
else:
    long_description = ''''''

setup(
    name='MolVS',
    version=molvs.__version__,
    author=molvs.__author__,
    author_email=molvs.__email__,
    license=molvs.__license__,
    url='https://github.com/mcs07/MolVS',
    packages=['molvs'],
    description='',
    long_description=long_description,
    keywords='chemistry cheminformatics rdkit',
    zip_safe=False,
    test_suite='nose.collector',
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
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
