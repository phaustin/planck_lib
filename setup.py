#!/usr/bin/env python
# https://docs.python.org/3/distutils/setupscript.html
# https://github.com/pypa/sampleproject/blob/main/setup.py
from setuptools import setup, find_packages

setup(
    name='planck_lib',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    classifiers=[
        'License :: OSI Approved :: BSD License'
    ],
    keywords='eoas ubc ocese',
    long_description="""e340 libraries""")
