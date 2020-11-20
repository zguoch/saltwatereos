from setuptools import find_packages, setup

from Cython.Build import cythonize
from distutils.extension import Extension

# build water
setup(name="pySWEOS",
    version="1.0.0",
    packages=find_packages(),
    author="Zhikui Guo",
    description="A python API of swEOS library.",
    long_description="swEOS is a C++ library designed to calculate phase relation and thermoproperties of NaCl-H2O (or say saltwater, seawater) system.",
    # long_description_content_type='text/markdown',
    url="https://gitlab.com/hydrothermal-openfoam/saltwatereos/-/tree/master/wrapper/python",
    ext_modules=cythonize(
    Extension(
    "H2O",
    language="c++",
    sources=["H2O.pyx","../../Library/src/H2O.C","../../Library/src/Fluid.C"],
    include_dirs=["../../Library/include"],
    )),
    install_requires=[
        'numpy;platform_system=="Darwin"',
        'numpy;platform_system=="Linux"'
    ]
    )