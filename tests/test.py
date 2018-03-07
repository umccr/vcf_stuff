import traceback
import os
import sys
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, expanduser
from datetime import datetime
import shutil
import subprocess

from ngs_utils.testing import BaseTestCase, check_call, vcf_ignore_lines, swap_output
from ngs_utils.utils import is_az, is_local, is_travis
from ngs_utils.file_utils import safe_mkdir, add_suffix, get_ungz_gz

from nose.plugins.attrib import attr


""" Prepare test data:
CHROM=16
REGIONS_CODE=16:0-100000
REGIONS_BED=16\\t0\\t100000
FASTA=/Users/vsaveliev/googledrive/bio/reference_data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
bcftools view -r $REGIONS_CODE ori/MB_100vs50-ensemble-annotated.vcf.gz -Oz -o test-ensemble.vcf.gz
bcftools view -r $REGIONS_CODE ori/MB_100vs50-mutect2-annotated.vcf.gz  -Oz -o test-mutect2.vcf.gz
bcftools view -r $REGIONS_CODE ori/MB_100vs50-vardict-annotated.vcf.gz  -Oz -o test-vardict.vcf.gz
bcftools view -r $REGIONS_CODE ori/MB_100vs50-strelka2-annotated.vcf.gz -Oz -o test-strelka2.vcf.gz
bcftools view -r $REGIONS_CODE ori/MB-benchmark.vcf.gz                  -Oz -o test-benchmark.vcf.gz
bedtools getfasta -fi $FASTA -bed <(echo "$REGIONS_BED") | sed "s/$REGIONS_CODE/$CHROM/" > test-GRCh37.fa
samtools faidx test-GRCh37.fa
"""

data_dir = join(dirname(__file__), BaseTestCase.data_dir)
input_ensemble_vcf = join(data_dir, 'test-ensemble.vcf.gz')
input_mutect2_vcf = join(data_dir, 'test-mutect2.vcf.gz')
input_vardict_vcf = join(data_dir, 'test-vardict.vcf.gz')
input_strelka2_vcf = join(data_dir, 'test-strelka2.vcf.gz')


@attr(kind='eval')
class TestEvalVcf(BaseTestCase):
    script = 'eval_vcf'
    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    ref_vcf = join(data_dir, 'test-benchmark.vcf.gz')
    genome_name = 'test-GRCh37'
    genome = join(data_dir, f'{genome_name}.fa')

    def setUp(self):
        BaseTestCase.setUp(self)

    def test_preset_ref(self):
        out_dir = join(TestEvalVcf.results_dir, 'preset_ref')
        cmdl = f'eval_vcf test-mb {input_ensemble_vcf} {input_vardict_vcf} ' \
               f'-g {TestEvalVcf.genome_name} -o {out_dir}'
        self._run_cmd(cmdl, [input_ensemble_vcf, input_vardict_vcf], out_dir)
        self._check_file_throws(join(out_dir, 'report.tsv'), ignore_matching_lines=vcf_ignore_lines)

    def test_custom_ref(self):
        out_dir = join(TestEvalVcf.results_dir, 'custom_ref')
        cmdl = f'eval_vcf {TestEvalVcf.ref_vcf} {input_ensemble_vcf} {input_vardict_vcf} ' \
               f'-g {TestEvalVcf.genome} -o {out_dir}'
        self._run_cmd(cmdl, [TestEvalVcf.ref_vcf, input_ensemble_vcf, input_vardict_vcf], out_dir)
        self._check_file_throws(join(out_dir, 'report.tsv'), ignore_matching_lines=vcf_ignore_lines)


@attr(kind='pon')
class TestPonAnno(BaseTestCase):
    script = 'pon_anno'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def test_pon_anno(self):
        out_vcf = join(TestPonAnno.results_dir, basename(add_suffix(input_vardict_vcf, 'pon')))
        cmdl = f'pon_anno {input_vardict_vcf} -o {out_vcf} -h 1'
        self._run_cmd(cmdl, input_vardict_vcf, out_vcf)
        self._check_file_throws(out_vcf, ignore_matching_lines=vcf_ignore_lines)

@attr(kind='pon')
class TestPoNPipeline(BaseTestCase):
    script = 'pon_pipeline'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def test_pon_pipeline(self):
        cmdl = f'pon_pipeline {input_strelka2_vcf} {input_vardict_vcf}' \
               f' -o {TestPoNPipeline.results_dir} -h1,2'
        self._run_cmd(cmdl, [input_strelka2_vcf, input_vardict_vcf], TestPoNPipeline.results_dir)


@attr(kind='norm')
class TestNormVcf(BaseTestCase):
    script = 'norm_vcf'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def _run_norm_vcf(self, input_vcf=None):
        out_vcf = join(TestNormVcf.results_dir, basename(add_suffix(input_vcf, 'norm')))
        cmdl = f'norm_vcf {input_vcf} -o {out_vcf} -g test-GRCh37'
        self._run_cmd(cmdl, [input_vcf], out_vcf)
        self._check_file_throws(out_vcf, ignore_matching_lines=vcf_ignore_lines)

    def test_norm_vcf_vardict(self):
        self._run_norm_vcf(input_vardict_vcf)

    def test_norm_vcf_strelka2(self):
        self._run_norm_vcf(input_strelka2_vcf)

    def test_norm_vcf_mutect(self):
        self._run_norm_vcf(input_mutect2_vcf)

@attr(kind='norm')
class TestPcgrPrep(BaseTestCase):
    script = 'pcgr_prep'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def _run_pcgr_prep(self, input_vcf=None):
        ungz, _ = get_ungz_gz(input_vcf)
        out_vcf = join(TestPcgrPrep.results_dir, basename(add_suffix(ungz, 'pcgr')))
        cmdl = f'pcgr_prep {input_vcf} > {out_vcf}'
        self._run_cmd(cmdl, [input_vcf], out_vcf)
        self._check_file_throws(out_vcf, ignore_matching_lines=vcf_ignore_lines)

    def test_pcgr_prep_vardict(self):
        self._run_pcgr_prep(input_vardict_vcf)

    def test_pcgr_prep_strelka2(self):
        self._run_pcgr_prep(input_strelka2_vcf)

    def test_pcgr_prep_mutect(self):
        self._run_pcgr_prep(input_mutect2_vcf)
