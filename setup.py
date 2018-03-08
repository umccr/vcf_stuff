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
    entry_points={
        'console_scripts': [
            'norm_vcf     = vcf_stuff.vcf_normalisation.norm_vcf:main',
            'pon_anno     = vcf_stuff.panel_of_normals.annotate:main',
            'pon_pipeline = vcf_stuff.panel_of_normals.pipeline:main',
            'pcgr_prep    = vcf_stuff.vcf_normalisation.pcgr_prep:main',
            'eval_vcf     = vcf_stuff.vcf_evaluation.eval_vcf:main',
        ],
    },
    include_package_data=True,
)
