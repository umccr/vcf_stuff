from os.path import isfile, join
from collections import OrderedDict

from ngs_utils.file_utils import verify_file
from ngs_utils.utils import OrderedDefaultDict


def iter_lines(fpath):
    if fpath:
        with open(fpath) as f:
            for l in f:
                l = l.replace('\n', '')
                if not l or l.startswith('#'):
                    continue
                yield l

def _read_list(reason, fpath):
    gene_d = {}
    fpath = verify_file(fpath, description=reason + ' blacklist genes file', is_critical=True)
    for l in iter_lines(fpath):
        fs = l.split('\t')
        gene_name = l.split('\t')[0]
        meta_info = l.split('\t')[1] if len(fs) == 2 else ''
        gene_d[gene_name] = meta_info
    return gene_d

def parse_gene_blacklists(genes_conf_section, incidentalome_dir):
    _d = OrderedDict()
    if 'published' in genes_conf_section:
        _d['Freq mut gene in HGMD'] = 'published/flags_in_hgmd.txt'
        _d['Freq mut gene in OMIM'] = 'published/flags_in_omim.txt'
        _d['Freq mut gene'] = 'published/flags.txt'
        _d['Incidentalome gene'] = 'published/incidentalome.txt'
        _d['MutSigCV gene'] = 'published/mutsigcv.txt'
    if 'low_complexity' in genes_conf_section:
        _d['Low complexity gene'] = 'low_complexity/low_complexity_entire_gene.txt'
    if 'repetitive_single_exome' in genes_conf_section:
        _d['Repetitive single exon gene'] = 'low_complexity/repetitive_single_exon_gene.txt'
    if 'abnormal_gc' in genes_conf_section:
        _d['Low GC gene'] = 'low_complexity/low_gc.txt'
        _d['High GC gene'] = 'low_complexity/high_gc.txt'
    if 'too_many_cosmic_mutations' in genes_conf_section:
        _d['Gene with too many COSMIC mutations'] = 'low_complexity/too_many_cosmic_mutations.txt'

    d = OrderedDefaultDict(dict)
    for reason, fn in _d.items():
        d[reason] = _read_list(reason, join(incidentalome_dir, fn))
    return d

def all_blacklisted_genes(genes_conf_section, incidentalome_dir):
    gene_blacklists_by_reason = parse_gene_blacklists(genes_conf_section, incidentalome_dir)
    return set(g for blacklist in gene_blacklists_by_reason.values() for g in blacklist.keys())
    

def parse_genes_list(fpath):
    genes = []
    if fpath and verify_file(fpath):
        genes = [line.strip() for line in open(fpath)]
    return genes

def check_gene_in_a_blacklist(gene_name, blacklist, aa_pos=None):
    if gene_name in blacklist:
        meta_info = blacklist[gene_name]
        if meta_info == '':
            return True
        elif aa_pos is not None:
            fs = meta_info.split(':')  # regions in form of :232, 553:, 42:111
            if fs[0] and aa_pos >= int(fs[0]):
                return True
            if fs[1] and aa_pos < int(fs[1]):
                return True
    return False

import yaml
from os.path import dirname, join, abspath, basename


def get_anno_config():
    with open(join(abspath(dirname(__file__)), 'anno_config.yaml')) as f:
        anno_cfg = yaml.load(f)
        return anno_cfg


def get_vcfanno_toml(work_dir, genome_cfg):
    toml_tmpl = join(abspath(dirname(__file__)), 'vcfanno.toml')

    with open(toml_tmpl) as f:
        tmpl_content = f.read()

    for dbname, dbpath in genome_cfg.items():
        tmpl_content = tmpl_content.replace('{{ ' + dbname + ' }}', dbpath)

    toml_final = join(work_dir, basename(toml_tmpl))
    with open(toml_final, 'w') as f:
        f.write(tmpl_content)

    return toml_final