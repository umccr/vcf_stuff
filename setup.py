#!/usr/bin/env python
import sys
import os
from os.path import join, isfile, abspath, dirname
from setuptools import setup

version = '0.1'
name = 'vcf_stuff'

setup(
    name=name,
    version=version,
    author='Vlad Saveliev',
    description='Evaluating, filtering, comparing, and visualising VCF',
    keywords='bioinformatics',
    license='GPLv3',
    packages=[
        name,
    ],
    scripts=[
        'scripts/panel_of_normals',
        'scripts/anno_pon',
        'scripts/normalise_vcf',
        'scripts/pcgr_prep',
        'scripts/eval_vcf',
    ],
    include_package_data=True,
)
