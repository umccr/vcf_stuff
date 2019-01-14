from __future__ import division
import os
from os.path import abspath, expanduser, join, dirname, pardir
from yaml import load as load_yaml

from ngs_utils.config import load_yaml_config, fill_dict_from_defaults
from ngs_utils.file_utils import verify_file, verify_module, adjust_path
from ngs_utils.logger import info, critical, warn, err


configs_dirpath = join(dirname(abspath(__file__)), 'configs')
filt_info_defaults_yaml_fpath = join(configs_dirpath, 'RUNINFO_DEFAULTS.yaml')
verify_file(filt_info_defaults_yaml_fpath, is_critical=True)
filt_info_defaults = load_yaml(open(filt_info_defaults_yaml_fpath))
filt_cnf_fpaths = dict(
    exome    = join(configs_dirpath, 'run_info_ExomeSeq.yaml'),
    genome   = join(configs_dirpath, 'run_info_WGS.yaml'),
    panel    = join(configs_dirpath, 'run_info_DeepSeq.yaml'),
    rnaseq   = join(configs_dirpath, 'run_info_RNAseq.yaml'),
)
filt_cnf_fpaths['targeted'] = filt_cnf_fpaths['deep_seq'] = filt_cnf_fpaths['panel']
filt_cnf_fpaths['wgs'] = filt_cnf_fpaths['genome']

def get_filt_cfg(opts, target_type=None, vardict_min_freq=None, is_wgs=False):
    if not isinstance(opts, dict):
        opts = opts.__dict__
    if not target_type and 'target_type' in opts:
        target_type = opts['target_type']

    # pick the defaults yaml based on target type, or from input filt_cnf
    filt_cnf = load_filt_cfg(
        verify_file(opts['filt_cnf'], is_critical=True) if opts['filt_cnf'] else None,
        target_type, vardict_min_freq, is_wgs)
    if not filt_cnf:
        return None
    
    # replace min_freq and act_min_freq from cmdl opts
    if 'min_freq' in opts and opts['min_freq'] is not None:
        filt_cnf['min_freq'] = float(opts['min_freq'])
        filt_cnf['act_min_freq'] = filt_cnf['min_freq'] / 2
        del filt_cnf['filt_cnf_fpath']
    filt_cnf['opt_d'] = opts
    return filt_cnf


def reload_filt_cfg(filt_cfg, target_type=None, vardict_min_freq=None, is_wgs=False, config_dir=None):
    if not filt_cfg:
        return None
    assert 'opt_d' in filt_cfg
    if config_dir is not None:
        run_info = detect_run_info_in_config_dir(config_dir)
        if run_info:
            filt_cfg['opt_d']['filt_cnf'] = run_info
    return get_filt_cfg(filt_cfg['opt_d'], target_type, vardict_min_freq, is_wgs)


def save_filt_cnf(filt_cnf, out_dir):
    filt_cnf_fpath = join(out_dir, 'run_info.yaml')
    with open(filt_cnf_fpath, 'w') as f:
        import yaml
        yaml.dump({k: v for k, v in filt_cnf.items()
                   if k not in ['blacklist', 'opt_d', 'filt_cnf_fpath']},
                  f, default_flow_style=False)
    return filt_cnf_fpath


def load_filt_cfg(filt_cnf_fpath=None, target_type=None, vardict_min_freq=None, is_wgs=False):
    """
    Specify either target_type, or vardict_min_freq and is_wgs
    """
    if not filt_cnf_fpath:
        if not target_type:
            if vardict_min_freq is not None:
                if vardict_min_freq <= 0.005:
                    info('Filtering config: min_allele_fraction=' + str(vardict_min_freq) + ' which is less 0.005, '
                         'setting config for panel')
                    target_type = 'panel'
                elif is_wgs is None:  # coverage interval is not defined
                    warn('Coverage interval is not defined, skipping variant filtering')
                    return None
                elif is_wgs is True:
                    target_type = 'genome'
                    info('Filtering config: setting config for genome')
                else:
                    target_type = 'exome'
                    info('Filtering config: min_allele_fraction=' + str(vardict_min_freq) + ' which is higher than 0.005, '
                         'setting config for exome')
            else:
                target_type = 'exome'
                info('Neither min freq not filt config was provided, using settings for exome')
        assert target_type in filt_cnf_fpaths, \
            'filt_cnf_fpath=' + str(filt_cnf_fpath) + '; ' + str(target_type) + ' not in ' + str(filt_cnf_fpaths.keys())
        filt_cnf_fpath = filt_cnf_fpaths[target_type]
    d = load_yaml_config(filt_cnf_fpath)
    if d.get('variant_filtering') and isinstance(d.get('variant_filtering'), dict):
        d = d.get('variant_filtering', dict())
    d = fill_dict_from_defaults(d, filt_info_defaults)
    d['filt_cnf_fpath'] = filt_cnf_fpath
    return d


def get_dbsnp_multi_mafs(genome_cfg):
    if 'dbsnp_multi_mafs' not in genome_cfg:
        warn('Warning: dbsnp_multi_mafs not provided in the system configuration file for the genome.')
        return None
    return verify_file(genome_cfg['dbsnp_multi_mafs'], is_critical=True,
                                   description='dbSNP multi mafs file in system configuration file')


def detect_run_info_in_config_dir(config_dir):
    run_info_fpaths_in_config = [
        abspath(join(config_dir, fname))
        for fname in os.listdir(config_dir)
        if fname.startswith('run_info') and fname.endswith('.yaml')]

    if len(run_info_fpaths_in_config) > 1:
        critical('More than one YAML file containing run_info in name found in the config '
                 'directory ' + config_dir + ': ' + ' '.join(run_info_fpaths_in_config))

    if len(run_info_fpaths_in_config) == 0:
        return None
        
    run_cnf = verify_file(run_info_fpaths_in_config[0], is_critical=True)
    info('Using run configuration from the config directory ' + run_cnf)
    return run_cnf

