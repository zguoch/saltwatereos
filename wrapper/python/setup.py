from setuptools import setup

from Cython.Build import cythonize
from distutils.extension import Extension

# build water
setup(ext_modules=cythonize(
    Extension(
    "H2O",
    language="c++",
    sources=["H2O.pyx","../../Library/src/H2O.C","../../Library/src/Fluid.C"],
    include_dir=["../../Library/include"],
    )))