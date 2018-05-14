#!/usr/bin/env python
import sys
import click
import os
import tempfile
from ngs_utils.call_process import run_simple
from python_utils.hpc import get_loc, find_loc, get_ref_file
import locale

from vcf_stuff.panel_of_normals import get_snps_toml_path, get_indels_toml_path

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
def main(vcf, genome, output_file=None, filter_hits=None):
    """
    Annotate records in the VCF file `vcf` with `INFO/PoN_CNT, optionally flag INFO/PoN_CNT>=filter_hits with FILTER=PoN`
    """

    fixed_snps_toml_f = tempfile.NamedTemporaryFile(delete=False)
    normals_dir = get_ref_file(genome, key='panel_of_normals_dir')
    run_simple(f'sed s#file=\\\"#file=\\\"{normals_dir}/# {get_snps_toml_path()} > {fixed_snps_toml_f.name}')

    run_simple(f'bcftools view -H {os.path.join(normals_dir, "panel_of_normals.snps.vcf.gz")}')

    fixed_indels_toml_f = tempfile.NamedTemporaryFile(delete=False)
    normals_dir = get_ref_file(genome, key='panel_of_normals_dir')
    run_simple(f'sed s#file=\\\"#file=\\\"{normals_dir}/# {get_indels_toml_path()} > {fixed_indels_toml_f.name}')

    run_simple(f'cat {fixed_indels_toml_f.name}')

    cmd = (f'vcfanno {fixed_snps_toml_f.name} {vcf}' +
        f' | vcfanno -permissive-overlap {fixed_indels_toml_f.name} /dev/stdin | bgzip -c')
    if filter_hits:
        cmd += f' | bcftools filter -e "INFO/PoN_CNT>={filter_hits}" --soft-filter PoN --mode + -Oz'

    if output_file:
        run_simple(cmd + f' > {output_file}')
        sys.stderr.write(f'Saved to {output_file}\n')
    else:
        run_simple(cmd)

    run_simple(f'bcftools view -H {output_file}')

    fixed_snps_toml_f.close()
    fixed_indels_toml_f.close()
    os.unlink(fixed_snps_toml_f.name)
    os.unlink(fixed_indels_toml_f.name)


if __name__ == '__main__':
    main()
