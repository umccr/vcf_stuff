from os.path import join, basename, splitext, dirname, abspath


def package_path():
    return dirname(abspath(__file__))


def get_gnomad_lua():
    return join(package_path(), 'anno_gnomad.lua')
