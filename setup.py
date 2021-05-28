#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 09:03:11 2020
@author: jguterl, bgarcia
"""

# Always prefer setuptools over distutils
from setuptools import setup,find_packages

import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='INGRID',
    version_format='{tag}',
    setup_requires=['setuptools-git-version'],
    description='Tokamak edge plasma grid generator',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='http://github.com/LLNL/Ingrid',
    author='B. Garcia, M. Umansky, J. Guterl',  # Optional

    # This should be a valid email address corresponding to the author listed
    # above.
    author_email='bgarci26@ucsc.edu, umansky1@llnl.gov',
    classifiers=[
	'Development Status :: 5 - Production/Stable',
        #'Intended Audience :: Any UEDGE users',
        'Topic :: Software Development',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
	'Programming Language :: Python :: 3.9',
    ],
    keywords='grid, mesh, generator, tokamak, edge, plasma, efit, uedge',
    packages=find_packages(),
    python_requires='>=3.5, <4',
    install_requires=[
        'm2r2',
        'scipy >= 1.3.1',
        'numpy',
        'pyyaml',
        'matplotlib',
        'tk',
        'sympy'
    ]
)

