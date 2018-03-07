#!/usr/bin/env python

import click
import os
import tempfile
from ngs_utils.call_process import run_simple
from python_utils.hpc import get_loc, find_loc, get_ref_file
import locale

from vcf_stuff.panel_of_normals import get_toml_path

try:
    if 'UTF-8' not in locale.getlocale(locale.LC_ALL):
        locale.setlocale(locale.LC_ALL, 'en_AU.UTF-8')
except TypeError:
    pass


@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('-g', 'genome', default='GRCh37')
@click.option('-o', 'output_file', type=click.Path())
@click.option('-h', 'filter_hits', type=click.INT, default=None)
@click.option('--check-allele', is_flag=True, default=False)
def main(vcf, genome, output_file=None, filter_hits=None, check_allele=False):
    """
    Annotate records in the VCF file `vcf` with `INFO/PoN_CNT`
    """
    fixed_toml_f = tempfile.NamedTemporaryFile(delete=False)

    normals_dir = get_ref_file(genome, key='panel_of_normals_dir')
    run_simple(f'sed s#file=\\\"#file=\\\"{normals_dir}/# {get_toml_path()} > {fixed_toml_f.name}')

    cmd = f'vcfanno{" " if check_allele else " -permissive-overlap"} {fixed_toml_f.name} {vcf} | bgzip -c'
    if filter_hits:
        cmd += f' | bcftools filter -e "INFO/PoN_CNT>={filter_hits}" --soft-filter PoN --mode + -Oz'

    if output_file:
        run_simple(cmd + f' > {output_file}')
        print('Saved to ' + output_file)
    else:
        run_simple(cmd)

    fixed_toml_f.close()
    os.unlink(fixed_toml_f.name)


if __name__ == '__main__':
    main()
