#!/usr/bin/env python
from os.path import join
import sys
import click
import os
import tempfile
from ngs_utils.call_process import run_simple
from python_utils.hpc import get_ref_file
from vcf_stuff.panel_of_normals import get_snps_toml_path, get_indels_toml_path
from ngs_utils.file_utils import verify_file, verify_dir
from ngs_utils.utils import set_locale; set_locale()

@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('-o', '--output-file', type=click.Path())
@click.option('-h', '--filter-hits', type=click.INT, default=None)
@click.option('-g', '--genome', default='GRCh37')
@click.option('--panel-of-normals-dir', '--pon-dir', default=None)
def main(vcf, output_file=None, filter_hits=None, genome=None, pon_dir=None):
    """
    Annotate records in the VCF file `vcf` with `INFO/PoN_CNT, optionally flag INFO/PoN_CNT>=filter_hits with FILTER=PoN`
    """

    fixed_snps_toml_f = tempfile.NamedTemporaryFile(delete=False)

    if pon_dir:
        verify_file(join(pon_dir, 'panel_of_normals.snps.vcf.gz'), is_critical=True, description='Panel of normals SNPs file in user provided folder')
        verify_file(join(pon_dir, 'panel_of_normals.snps.vcf.gz.tbi'), is_critical=True, description='Please index panel of normal files with tabix')
        verify_file(join(pon_dir, 'panel_of_normals.indels.vcf.gz'), is_critical=True, description='Panel of normals indels file in user provided folder')
        verify_file(join(pon_dir, 'panel_of_normals.indels.vcf.gz.tbi'), is_critical=True, description='Please index panel of normal files with tabix')
        pon_dir = verify_dir(pon_dir, is_critical=True)
    else:
        pon_dir = get_ref_file(genome, key='panel_of_normals_dir')

    run_simple(f'sed s#file=\\\"#file=\\\"{pon_dir}/# {get_snps_toml_path()} > {fixed_snps_toml_f.name}')

    # run_simple(f'bcftools view -H {os.path.join(normals_dir, "panel_of_normals.snps.vcf.gz")}')

    fixed_indels_toml_f = tempfile.NamedTemporaryFile(delete=False)
    run_simple(f'sed s#file=\\\"#file=\\\"{pon_dir}/# {get_indels_toml_path()} > {fixed_indels_toml_f.name}')

    # run_simple(f'cat {fixed_indels_toml_f.name}')

    cmd = (f'vcfanno {fixed_snps_toml_f.name} {vcf}' +
        f' | vcfanno -permissive-overlap {fixed_indels_toml_f.name} /dev/stdin | bgzip -c')
    if filter_hits:
        cmd += f' | bcftools filter -e "INFO/PoN_CNT>={filter_hits}" --soft-filter PoN --mode +'

    if output_file:
        run_simple(cmd + f' -Oz -o {output_file}')
        sys.stderr.write(f'Saved to {output_file}\n')
    else:
        run_simple(cmd)

    # run_simple(f'bcftools view -H {output_file}')

    fixed_snps_toml_f.close()
    fixed_indels_toml_f.close()
    os.unlink(fixed_snps_toml_f.name)
    os.unlink(fixed_indels_toml_f.name)


if __name__ == '__main__':
    main()
