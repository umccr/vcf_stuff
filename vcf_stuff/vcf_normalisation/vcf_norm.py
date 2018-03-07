#!/usr/bin/env python
"""
Normalizes VCF: 
- splits multiallelic ALT
- splits biallelic MNP
- left-aligns indels
- fixes FORMAT and INFO fields
"""
import os
import click
from os.path import isfile

from python_utils.hpc import get_ref_file
from ngs_utils.call_process import run
from ngs_utils.file_utils import verify_file


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', 'output_file', type=click.Path())
@click.option('-g', '-f', 'genome', type=click.Path(), default='GRCh37', help='Path to genome fasta, or genome build name (if a known location)')
def main(input_file, output_file, genome=False):
    reference_fasta = get_ref_file(genome)
    cmd = make_normalise_cmd(input_file, output_file, reference_fasta)
    run(cmd)


def make_normalise_cmd(input_file, output_file, reference_fasta):
    return (
        f'bcftools norm -m \'-\' {input_file} -Ov -f {reference_fasta}'
        f' | vcfallelicprimitives -t DECOMPOSED --keep-geno --keep-info'
        f' | vcfstreamsort'
        f' | grep -v "##INFO=<ID=TYPE,Number=1"'
        f' | bgzip -c > {output_file}'
        f' && tabix -p vcf {output_file}')


if __name__ == '__main__':
    main()

