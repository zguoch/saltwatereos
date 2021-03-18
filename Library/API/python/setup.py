#!/usr/bin/env python
# coding: utf-8

from setuptools import setup, find_packages
version=1.0
try:
    with open('version.txt') as f:
        lines = f.readlines()
        version=lines[0]
except:
    print('Please provide version file in **/API/python folder')
    exit(0)
    pass
setup(
    name='pyswEOS',
    version=version,
    author='Zhikui Guo',
    author_email='zguo@geomar.de',
    url='https://www.scibyte.cn',
    python_requires='==3.8',
    description=u'Python API of swEOS library which is designed to calculate EOS and thermal dynamic properties of H2O-NaCl system, usually called salt water.',
    packages=find_packages(where='.', exclude=('docs_pyH2ONaCl'), include=('*',)), 
    install_requires=[
        'matplotlib',
        'numpy',
        'argparse',
        'colored'
        ],
    include_package_data=True,
    entry_points={
        'console_scripts': []
    },
    keywords = "NaCl H2O IAPWS Hydrothermal seawater"
)