from os.path import dirname, abspath, join, basename, isfile
from ngs_utils.file_utils import verify_file, splitext_plus


def package_path():
    return dirname(abspath(__file__))


def get_snps_toml_path():
    return verify_file(join(package_path(), 'vcfanno_snps.toml'))


def get_indels_toml_path():
    return verify_file(join(package_path(), 'vcfanno_indels.toml'))

