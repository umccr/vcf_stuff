#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup, find_packages
import releazit

import vcf_stuff
pkg = vcf_stuff.__name__

version = releazit.get_version(pkg)

setup(
    name=pkg,
    version=version,
    author='Vlad Saveliev',
    author_email='vladislav.sav@gmail.com',
    description='Evaluating, filtering, comparing, and visualising variant calls',
    keywords='bioinformatics',
    url='https://github.com/umccr/vcf-stuff',
    license='GPLv3',
    packages=find_packages(),
    package_data={
        pkg: releazit.find_package_files('', pkg)
    },
    scripts=[join('scripts', fn) for fn in os.listdir('scripts')],
    include_package_data=True,
    zip_safe=False,
)
