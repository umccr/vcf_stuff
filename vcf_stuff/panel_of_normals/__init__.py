from os.path import dirname, abspath, join, basename, isfile
import click
import os
import sys
import yaml
import tempfile
import subprocess
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import verify_file, splitext_plus
from python_utils.hpc import get_loc, find_loc, get_ref_file
import locale

try:
    if 'UTF-8' not in locale.getlocale(locale.LC_ALL):
        locale.setlocale(locale.LC_ALL, 'en_AU.UTF-8')
except TypeError:
    pass


def package_path():
    return dirname(abspath(__file__))


def get_toml_path():
    return verify_file(join(package_path(), 'vcfanno.toml'))


@click.command()
@click.argument('vcf', type=click.Path(exists=True))
@click.option('-g', 'genome', default='GRCh37')
@click.option('-o', 'output_file', type=click.Path())
@click.option('-h', 'filter_hits', type=click.INT, default=None)
@click.option('--check-allele', is_flag=True, default=False)
def annotate(vcf, genome, output_file=None, filter_hits=None, check_allele=False):
    """
    Annotate records in the VCF file `vcf` with `INFO/PoN_CNT`
    """
    fixed_toml_f = tempfile.NamedTemporaryFile(delete=False)

    normals_dir = get_ref_file(genome, key='panel_of_normals_dir')
    toml = get_toml_path()
    run_simple(f'sed s#file=\\\"#file=\\\"{normals_dir}/# {toml} > {fixed_toml_f.name}')

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


@click.command()
@click.argument('vcfs', nargs=-1, type=click.Path(exists=True))
@click.option('-g', 'genome', default='GRCh37')
@click.option('-o', 'output_dir', type=click.Path())
@click.option('-j', 'jobs', type=click.INT, default=1)
@click.option('-h', 'hits_thresholds')
def pipeline(vcfs, genome, output_dir=None, jobs=1, hits_thresholds=None):
    """
    Filter all VCF files `vcfs`: remove records with `PoN_CNT` > `hits_thresholds`
    """
    config = {
        'samples': {splitext_plus(basename(v))[0]: abspath(v) for v in vcfs
                    if v.endswith('.vcf') or v.endswith('.vcf.gz')},
        'hits_thresholds': hits_thresholds.split(',') if hits_thresholds else [1, 2, 3],
        'genome': genome,
    }

    f = tempfile.NamedTemporaryFile(mode='wt', delete=False)
    yaml.dump(config, f)
    f.close()

    cmd = (f'snakemake ' +
           f'--snakefile {join(package_path(), "Snakefile")} ' +
           f'--printshellcmds ' +
          (f'--directory {output_dir} ' if output_dir else ' ') +
           f'--configfile {f.name} '
           f'--jobs {jobs} ')
    print(cmd)
    os.system(cmd)
