#!/usr/bin/env python
# coding: utf-8

from setuptools import setup, find_packages

setup(
    name='pyH2ONaCl',
    version='0.0.1',
    author='Zhikui Guo',
    author_email='zguo@geomar.de',
    url='https://www.scibyte.cn',
    description=u'Python API of H2ONaCl library which is designed to calculate EOS and thermal dynamic properties of H2O-NaCl system, usually called salt water.',
    packages=find_packages(where='.', exclude=('docs_pyH2ONaCl'), include=('*',)), 
    install_requires=[
        'matplotlib',
        'numpy'
        ],
    entry_points={
        'console_scripts': []
    },
    keywords = "NaCl H2O IAPWS Hydrothermal seawater"
)