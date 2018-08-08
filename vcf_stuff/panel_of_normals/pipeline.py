#!/usr/bin/env python

from os.path import dirname, abspath, join, basename, isfile
import click
import os
import sys
import yaml
import tempfile
from ngs_utils.file_utils import splitext_plus
from hpc_utils.hpc import get_ref_file
from vcf_stuff.panel_of_normals import package_path

from ngs_utils.utils import set_locale
set_locale()

@click.command()
@click.argument('vcfs', nargs=-1, type=click.Path(exists=True))
@click.option('-g', 'genome', default='GRCh37')
@click.option('-o', 'output_dir', type=click.Path())
@click.option('-j', 'jobs', type=click.INT, default=1)
@click.option('-h', 'hits_thresholds')
def main(vcfs, genome, output_dir=None, jobs=1, hits_thresholds=None):
    """
    Filter all VCF files `vcfs`: remove records with `PoN_CNT` > `hits_thresholds`
    """
    ref_fa = get_ref_file(genome)
    normals_dir = get_ref_file(genome, key='panel_of_normals_dir')

    config = {
        'samples': {splitext_plus(basename(v))[0]: abspath(v) for v in vcfs
                    if v.endswith('.vcf') or v.endswith('.vcf.gz')},
        'hits_thresholds': hits_thresholds.split(',') if hits_thresholds else [1, 2, 3],
        'ref_fa': ref_fa,
        'normals_dir': normals_dir,
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


if __name__ == '__main__':
    main()
