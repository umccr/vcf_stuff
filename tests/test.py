from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, abspath, expanduser
from nose.plugins.attrib import attr
from ngs_utils.testing import BaseTestCase, check_call, vcf_ignore_lines, swap_output
from ngs_utils.file_utils import safe_mkdir, add_suffix, get_ungz_gz



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

ssh spa
bcftools view -r $REGIONS_CODE /data/cephfs/punim0010/extras/panel_of_normals/panel_of_normals.snps.vcf.gz -Oz -o panel_of_normals.TEST.snps.vcf.gz
bcftools view -r $REGIONS_CODE /data/cephfs/punim0010/extras/panel_of_normals/panel_of_normals.indels.vcf.gz -Oz -o panel_of_normals.TEST.indels.vcf.gz
tabix -p vcf panel_of_normals.TEST.snps.vcf.gz
tabix -p vcf panel_of_normals.TEST.indels.vcf.gz

mkdir panel_of_normals
scp -r spa:/data/cephfs/punim0010/extras/panel_of_normals/panel_of_normals/panel_of_normals.TEST.indels.vcf.gz panel_of_normals
scp -r spa:/data/cephfs/punim0010/extras/panel_of_normals/panel_of_normals/panel_of_normals.TEST.indels.vcf.gz.tbi panel_of_normals
scp -r spa:/data/cephfs/punim0010/extras/panel_of_normals/panel_of_normals/panel_of_normals.TEST.snps.vcf.gz panel_of_normals
scp -r spa:/data/cephfs/punim0010/extras/panel_of_normals/panel_of_normals/panel_of_normals.TEST.snps.vcf.gz.tbi panel_of_normals

bedtools getfasta -fi $FASTA -bed <(echo "$REGIONS_BED") | sed "s/$REGIONS_CODE/$CHROM/" > test-GRCh37.fa
samtools faidx test-GRCh37.fa
"""


data_dir = join(dirname(__file__), BaseTestCase.data_dir)
input_ensemble_vcf = join(data_dir, 'test-ensemble.vcf.gz')
input_mutect2_vcf = join(data_dir, 'test-mutect2.vcf.gz')
input_vardict_vcf = join(data_dir, 'test-vardict.vcf.gz')
input_strelka2_vcf = join(data_dir, 'test-strelka2.vcf.gz')
genome_name = 'test-GRCh37'
ref_fa = join(data_dir, f'{genome_name}.fa')


@attr(kind='eval')
class TestEvalVcf(BaseTestCase):
    script = 'eval_vcf'
    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    ref_vcf = join(data_dir, 'test-benchmark.vcf.gz')

    def setUp(self):
        BaseTestCase.setUp(self)

    def test_preset_ref(self):
        out_dir = join(TestEvalVcf.results_dir, 'preset_ref')
        cmdl = f'eval_vcf test-mb {input_ensemble_vcf} {input_vardict_vcf} ' \
               f'-g {genome_name} -o {out_dir}'
        self._run_cmd(cmdl, [input_ensemble_vcf, input_vardict_vcf], out_dir)
        self._check_file_throws(join(out_dir, 'report.tsv'), ignore_matching_lines=vcf_ignore_lines)

    def test_custom_ref(self):
        out_dir = join(TestEvalVcf.results_dir, 'custom_ref')
        cmdl = f'eval_vcf {TestEvalVcf.ref_vcf} {input_ensemble_vcf} {input_vardict_vcf} ' \
               f'--ref-fasta {ref_fa} -o {out_dir} -j4'
        self._run_cmd(cmdl, [TestEvalVcf.ref_vcf, input_ensemble_vcf, input_vardict_vcf], out_dir)
        self._check_file_throws(join(out_dir, 'report.tsv'), ignore_matching_lines=vcf_ignore_lines)


# @attr(kind='eval_cnv')
# class TestEvalCnv(BaseTestCase):
#     script = 'eval_cnv'
#     data_dir = join(dirname(__file__), BaseTestCase.data_dir, 'cnv', 'hcc2218')
#     results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
#     gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)
#
#     def setUp(self):
#         BaseTestCase.setUp(self)
#         self.ref_cnv = join(TestEvalCnv.data_dir, 'HCC2218_truthset_cnv_bcbio.tsv')
#         self.input_cnvs = [join(TestEvalCnv.data_dir, fn) for fn in [
#             'HCC2218_cnvkit-call.cns',
#             'HCC2218_facets_cncf.tsv',
#             'HCC2218_manta.vcf',
#             'HCC2218_purple.cnv.tsv',
#         ]]
#
#     def test(self):
#         out_dir = join(TestEvalCnv.results_dir)
#         cmdl = f'eval_cnv {self.ref_cnv} {" ".join(self.input_cnvs)} -g GRCh37 -o {out_dir}'
#         self._run_cmd(cmdl, self.input_cnvs, out_dir)
#         self._check_file_throws(join(out_dir, 'report.tsv'))
#         self._check_file_throws(join(out_dir, 'table.tsv'))


pon_data_dir = join(data_dir, 'panel_of_normals')

@attr(kind='pon')
class TestPonAnno(BaseTestCase):
    script = 'pon_anno'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def test_pon_anno(self):
        out_vcf = join(TestPonAnno.results_dir, basename(add_suffix(input_vardict_vcf, 'pon')))
        cmdl = f'pon_anno {input_vardict_vcf} -o {out_vcf} -h 1 -g {genome_name} --pon-dir {pon_data_dir}'
        self._run_cmd(cmdl, input_vardict_vcf, out_vcf)
        self._check_file_throws(out_vcf, ignore_matching_lines=vcf_ignore_lines)

@attr(kind='pon')
class TestPoNPipeline(BaseTestCase):
    script = 'pon_pipeline'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def test_pon_pipeline(self):
        cmdl = f'pon_pipeline {input_strelka2_vcf} {input_vardict_vcf}' \
               f' -o {TestPoNPipeline.results_dir} -h1,2 -g {genome_name} --pon-dir {pon_data_dir}'
        self._run_cmd(cmdl, [input_strelka2_vcf, input_vardict_vcf], TestPoNPipeline.results_dir)
        self._check_file_throws(join(TestPoNPipeline.results_dir, 'pon_filter', 'test-strelka2-n1.vcf.gz'), ignore_matching_lines=vcf_ignore_lines)
        self._check_file_throws(join(TestPoNPipeline.results_dir, 'pon_filter', 'test-vardict-n2.vcf.gz'), ignore_matching_lines=vcf_ignore_lines)


@attr(kind='norm')
class TestNormVcf(BaseTestCase):
    script = 'norm_vcf'
    results_dir = join(dirname(__file__), BaseTestCase.results_dir, script)
    gold_standard_dir = join(dirname(__file__), BaseTestCase.gold_standard_dir, script)

    def _run_norm_vcf(self, input_vcf=None):
        out_vcf = join(TestNormVcf.results_dir, basename(add_suffix(input_vcf, 'norm')))
        cmdl = f'norm_vcf {input_vcf} -o {out_vcf} --ref-fasta {ref_fa}'
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
        self._check_file_throws(out_vcf,
            wrapper='bcftools query -f "%TUMOR_AF\\t%NORMAL_AF\\t%TUMOR_DP\\t%NORMAL_DP\\t%TUMOR_MQ\\n" | '
                    'awk \'{printf("%.2f %.2f %.2f %.2f %.2f\\n", $1, $2, $3, $4, $5) }\'')

    def test_pcgr_prep_vardict(self):
        self._run_pcgr_prep(input_vardict_vcf)

    def test_pcgr_prep_strelka2(self):
        self._run_pcgr_prep(input_strelka2_vcf)

    def test_pcgr_prep_mutect(self):
        self._run_pcgr_prep(input_mutect2_vcf)
