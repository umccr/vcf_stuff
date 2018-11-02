#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup, find_packages
from ngs_utils import setup_utils

name = 'vcf_stuff'

version = setup_utils.get_cur_version(name)

setup(
    name=name,
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='Evaluating, filtering, comparing, and visualising variant calls',
    keywords='bioinformatics',
    url='https://github.com/umccr/vcf-stuff',
    license='GPLv3',
    package_data={
        name: setup_utils.find_package_files('', name)
    },
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    scripts=[join('scripts', fn) for fn in os.listdir('scripts')],
)
