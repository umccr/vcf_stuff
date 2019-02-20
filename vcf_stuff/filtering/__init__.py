from os.path import join, basename, splitext, dirname, abspath


def package_path():
    return dirname(abspath(__file__))


def get_gnomad_lua():
    return join(package_path(), 'anno_gnomad.lua')


def add_cyvcf2_filter(rec, filt):
    filters = rec.FILTER.split(';') if rec.FILTER else []
    filters = set(filters)
    filters.add(filt)
    f = ';'.join(filters)
    rec.FILTER = str(f)
    return rec
