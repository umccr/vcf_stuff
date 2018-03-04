from os.path import isfile, join, dirname, abspath
from ngs_utils.file_utils import verify_file


here = dirname(abspath(__file__))


def package_path():
    return here


def get_toml_path():
    return verify_file(join(here, 'vcfanno.toml'))
