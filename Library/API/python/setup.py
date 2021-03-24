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
AUTHOR = "Zhikui Guo"
AUTHOR_EMAIL = "zguo@geomar.de"
MAINTAINER = AUTHOR
MAINTAINER_EMAIL = AUTHOR_EMAIL
URL = "https://www.sweos.info"
with open("readme.rst") as f:
    LONG_DESCRIPTION = "".join(f.readlines())
CLASSIFIERS = [
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    "Intended Audience :: Education",
    "Topic :: Scientific/Engineering",
    "Topic :: Software Development :: Libraries",
    "Programming Language :: Python :: 2.7",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9"
]
setup(
    name='pyswEOS',
    version=version,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=MAINTAINER_EMAIL,
    url='https://www.sweos.info',
    description=u'Python API of swEOS library which is designed to calculate EOS and thermal dynamic properties of H2O-NaCl system, usually called salt water.',
    long_description=LONG_DESCRIPTION,
    classifiers=CLASSIFIERS,
    packages=find_packages(where='.', exclude=('docs_pyH2ONaCl'), include=('*',)), 
    install_requires=['matplotlib','numpy','argparse','colored'],
    include_package_data=True,
    entry_points={
        'console_scripts': []
    },
    keywords = "NaCl H2O IAPWS Hydrothermal seawater"
)