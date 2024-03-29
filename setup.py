#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup, find_packages
import versionpy

pkg = 'vcf_stuff'

setup(
    name=pkg,
    version='0.6.4', # versioned by bump2version
    author='Vlad Savelyev',
    author_email='vladislav.sav@gmail.com',
    description='Evaluating, filtering, comparing, and visualising variant calls',
    keywords='bioinformatics',
    url='https://github.com/umccr/vcf_stuff',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        pkg: versionpy.find_package_files('', pkg)
    },
    scripts=[join('scripts', fn) for fn in os.listdir('scripts')],
    include_package_data=True,
    zip_safe=False,
)
