#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup
from ngs_utils import setup_utils

name = 'vcf_stuff'

version = setup_utils.get_cur_version(name)

setup(
    name=name,
    version=version,
    author='Vlad Saveliev',
    description='Evaluating, filtering, comparing, and visualising variant calls',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        name,
    ],
    scripts=[join('scripts', fn) for fn in os.listdir('scripts')],
    include_package_data=True,
)
