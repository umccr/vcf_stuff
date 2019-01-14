from os.path import dirname, join, abspath

_p = abspath(dirname(__file__))


def suppressors():
    return join(_p, 'suppressors.txt')

def oncogenes():
    return join(_p, 'oncogenes.txt')

def incidentalome_dir():
    return join(_p, 'incidentalome')

#############
# Filtering
#
def _get_filt_file(fname, genome=None):
    short_genome = genome.split('-')[0] if genome else None
    return join(_p, short_genome or 'common', fname)

def common_snp(genome):
    return _get_filt_file('filter_common_snp.txt', genome)

def common_art(genome):
    return _get_filt_file('filter_common_artifacts.txt', genome)

def actionable(genome):
    return _get_filt_file('actionable.txt', genome)

def compendia(genome):
    return _get_filt_file('Compendia.MS7.Hotspot.txt', genome)

def splice(genome):
    return _get_filt_file('SPLICE.txt', genome)

def snpeffect_export_polymorphic():
    return _get_filt_file('snpeffect_export_Polymorphic.txt')

def actionable_hotspot():
    return _get_filt_file('actionable_hotspot.txt')

def ruledir():
    return _get_filt_file('rules')

def specific_mutations():
    return _get_filt_file('specific_mutations.tsv')

def last_critical_aa():
    return _get_filt_file('last_critical_aa.txt')

def ngs_reports_comments():
    return _get_filt_file('ngs_reports_comments.tsv')