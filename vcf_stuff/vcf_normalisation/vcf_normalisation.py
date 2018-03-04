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
from ngs_utils.file_utils import verify_file
from os.path import isfile
from umccrise.utils import get_loc

@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', 'output_file', type=click.Path())
@click.option('-f', 'reference_fasta', type=click.Path())
def main(input_file, output_file, reference_fasta=False):
    if not isfile(reference_fasta):
        reference_fasta = os.path.join(get_loc().hsapiens, reference_fasta)
    verify_file(reference_fasta, is_critical=True)
    cmd = make_normalise_cmd(input_file, output_file, reference_fasta)
    print(cmd)
    os.system(cmd)


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

