# -*- coding: utf-8 -*-
"""
Setup Script for LD_Direction
"""
import os
import codecs

from setuptools import setup

from Cython.Build import cythonize

###############################################################################
#                     Build the things we need for setup                      #
###############################################################################

# Get the long description from the README file
here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

###############################################################################
#                                Setup Options                                #
###############################################################################

setup(
    name='LD_Direction',
    version='0.1.0a',
    description=('Lookup LD between any two lists of SNPs'),
    long_description=long_description,
    url='https://github.com/MikeDacre/dbSNP',
    author='Michael Dacre',
    author_email='mike.dacre@gmail.com',
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Operating System :: Linux',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    keywords='likage-disequilibrium LD variation',

    install_requires=['tqdm', 'dbSNP', 'pandas'],
    packages=['LD_Direction'],
    ext_modules = cythonize('LD_Direction/distance_filtering.pyx'),
)
