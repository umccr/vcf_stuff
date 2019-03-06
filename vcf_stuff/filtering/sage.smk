from os.path import isfile, join, basename, splitext
from hpc_utils import hpc
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.vcf_utils import get_sample_ids, get_sample_names
from vcf_stuff import iter_vcf
import cyvcf2
from vcf_stuff.filtering import add_cyvcf2_filter, package_path


localrules: sage


SAMPLE = config['sample']
GENOME = config['genome']
TUMOR_BAM = config['tumor_bam']
NORMAL_BAM = config['normal_bam']
EXISTING_VCF = config.get('existing_vcf', None)
OUTPUT_VCF = config['output_vcf']
assert EXISTING_VCF.endswith('.vcf.gz'), EXISTING_VCF
assert OUTPUT_VCF.endswith('.vcf.gz'), OUTPUT_VCF
if EXISTING_VCF:
    TUMOR_NAME, NORMAL_NAME = get_sample_names(EXISTING_VCF)
else:
    TUMOR_NAME = splitext(basename(TUMOR_BAM))
    NORMAL_NAME = splitext(basename(NORMAL_BAM))

if config.get('genomes_dir'):
    hpc.genomes_dir = config.get('genomes_dir')


rule all:
    input:
        vcf = OUTPUT_VCF,
        tbi = OUTPUT_VCF + '.tbi'


rule run_sage:
    input:
        tumor_bam    = TUMOR_BAM,
        normal_bam   = NORMAL_BAM,
        coding_bed   = hpc.get_ref_file(GENOME, key='coding_regions'),
        ref_fa       = hpc.get_ref_file(GENOME, key='fa'),
        hotspots_vcf = hpc.get_ref_file(GENOME, key='hotspots'),
    output:
        sage_vcf = f'work/call/{SAMPLE}-sage.vcf.gz',
        sage_tbi = f'work/call/{SAMPLE}-sage.vcf.gz.tbi',
    params:
        jar = join(package_path(), 'sage-1.0-jar-with-dependencies.jar'),
        normal_sname = NORMAL_NAME,
        tumor_sname  = TUMOR_NAME,
        xms = 2000,
        xmx = 19000,
    resources:
        mem_mb = 20000
    group: "sage"
    shell:
        'java -Xms{params.xms}m -Xmx{params.xmx}m -cp {params.jar} com.hartwig.hmftools.sage.SageHotspotApplication '
        '-tumor {params.tumor_sname} -tumor_bam {input.tumor_bam} '
        '-reference {params.normal_sname} -reference_bam {input.normal_bam} '
        '-known_hotspots <(bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" {input.hotspots_vcf}) '
        '-coding_regions {input.coding_bed} '
        '-ref_genome {input.ref_fa} '
        '-out {output.sage_vcf} '
        '&& tabix -f -p vcf {output.sage_vcf}'


rule sage_rename_anno:
    input:
        sage_vcf = rules.run_sage.output.sage_vcf,
    output:
        sage_vcf = f'work/rename_anno/{SAMPLE}-sage.vcf.gz',
        sage_tbi = f'work/rename_anno/{SAMPLE}-sage.vcf.gz.tbi',
    group: "sage"
    shell: """
bcftools view {input.sage_vcf} | \
sed 's/HOTSPOT/SAGE_HOTSPOT/g' | \
sed 's/Hotspot Type: known, inframe/SAGE Hotspot Type: known, inframe'/g | \
bgzip -c > {output.sage_vcf} && tabix -p vcf {output.sage_vcf}
"""


rule sage_reorder_samples:
    input:
        sage_vcf = rules.sage_rename_anno.output.sage_vcf,
        sage_tbi = rules.sage_rename_anno.output.sage_tbi,
        vcf      = EXISTING_VCF,
    output:
        sage_vcf = f'work/sage_reorder_samples/{SAMPLE}-sage.vcf.gz',
        sage_tbi = f'work/sage_reorder_samples/{SAMPLE}-sage.vcf.gz.tbi',
    group: "sage"
    run:
        tumor_index, normal_index = get_sample_ids(input.vcf)
        assert sorted([tumor_index, normal_index]) == [0, 1]
        sample_in_order = [None, None]
        tumor_name, normal_name = get_sample_names(input.vcf)
        sample_in_order[tumor_index] = tumor_name
        sample_in_order[normal_index] = normal_name
        shell(f'bcftools view -s {",".join(sample_in_order)} {input.sage_vcf} -Oz -o {output.sage_vcf} '
              f'&& tabix -p vcf {output.sage_vcf}')


rule sage_pass:
    input:
        sage_vcf = rules.sage_reorder_samples.output.sage_vcf if EXISTING_VCF else rules.sage_rename_anno.output.sage_vcf,
        sage_tbi = rules.sage_reorder_samples.output.sage_tbi if EXISTING_VCF else rules.sage_rename_anno.output.sage_tbi,
    output:
        sage_vcf = f'work/sage_pass/{SAMPLE}-sage.vcf.gz',
        sage_tbi = f'work/sage_pass/{SAMPLE}-sage.vcf.gz.tbi',
    group: "sage"
    shell:
        'bcftools view -f.,PASS {input.sage_vcf} -Oz -o {output.sage_vcf} '
        '&& tabix -p vcf {output.sage_vcf}'


rule sage_pass_novel:
    input:
        sage_vcf = rules.sage_pass.output.sage_vcf,
        sage_tbi = rules.sage_pass.output.sage_tbi,
        vcf      = EXISTING_VCF,
    output:
        sage_vcf = f'work/sage_pass_novel/{SAMPLE}-sage.vcf.gz',
        sage_tbi = f'work/sage_pass_novel/{SAMPLE}-sage.vcf.gz.tbi',
    group: "sage"
    shell:
        'bcftools isec {input.sage_vcf} {input.vcf} -C -w1 -Oz -o {output.sage_vcf} '
        '&& tabix -p vcf {output.sage_vcf}'


rule add_novel_sage_calls:
    input:
        vcf      = EXISTING_VCF,
        sage_vcf = rules.sage_pass_novel.output.sage_vcf,
        sage_tbi = rules.sage_pass_novel.output.sage_tbi,
    output:
        vcf = f'work/add_novel_sage_calls/{SAMPLE}.vcf.gz',
        tbi = f'work/add_novel_sage_calls/{SAMPLE}.vcf.gz.tbi',
    group: "sage"
    run:
         shell('bcftools concat -a {input.vcf} {input.sage_vcf} -Oz -o {output.vcf} '
               '&& tabix -p vcf {output.vcf}')
         assert len(cyvcf2.VCF(output.vcf).samples) == 2

rule sort_saged:
    input:
        vcf = rules.add_novel_sage_calls.output.vcf,
    output:
        vcf = f'work/sort_saged/{SAMPLE}.vcf.gz',
        tbi = f'work/sort_saged/{SAMPLE}.vcf.gz.tbi',
    group: "sage"
    shell:
        '(bcftools view -h {input.vcf} ; bcftools view -H -f.,PASS {input.vcf} | sort -k1,1V -k2,2n) | '
        'bgzip -c > {output.vcf} && tabix -f -p vcf {output.vcf}'

rule annotate_from_sage:
    input:
        vcf = rules.sort_saged.output.vcf,
        tbi = rules.sort_saged.output.tbi,
        sage_vcf = rules.sage_reorder_samples.output.sage_vcf,
        sage_tbi = rules.sage_reorder_samples.output.sage_tbi,
    output:
        vcf = f'work/annotate_from_sage/{SAMPLE}.vcf.gz',
        tbi = f'work/annotate_from_sage/{SAMPLE}.vcf.gz.tbi',
    group: "sage"
    run:
        sage_calls = dict()
        for rec in cyvcf2.VCF(input.sage_vcf):
            key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
            sage_calls[key] = rec

        def proc_hdr(vcf):
            vcf.add_filter_to_header({'ID': 'SAGE_lowconf', 'Description': 'SAGE assigned low confidence to this call'})

        def proc_rec(rec, vcf, tumor_index, normal_index):
            key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
            sage_call = sage_calls.get(key)
            if sage_call is not None:
                # pcgr_prep should handle the existing SAGE fields, so we just need to
                # figure out how to populate FORMAT fields from cyvcf2
                if not sage_call.FILTER or sage_call.FILTER == 'PASS':
                    assert sage_call.INFO.get('SAGE_HOTSPOT') is not None, sage_call
                    rec.INFO['SAGE_HOTSPOT'] = sage_call.INFO['SAGE_HOTSPOT']
                    rec.FILTER = 'PASS'
                else:
                    add_cyvcf2_filter(rec, 'SAGE_lowconf')
                rec.set_format('DP', sage_call.format('DP'))
                rec.set_format('AD', sage_call.format('AD'))
            return rec

        tumor_index, normal_index = get_sample_ids(input.vcf)
        iter_vcf(input.vcf, output.vcf, proc_rec, proc_hdr=proc_hdr, tumor_index=tumor_index, normal_index=normal_index)

rule sage:
    input:
        rules.annotate_from_sage.output.vcf if EXISTING_VCF else rules.sage_pass_novel.output.sage_vcf,
    output:
        vcf = OUTPUT_VCF,
        tbi = OUTPUT_VCF + '.tbi'
    shell:
        'cp {input} {output.vcf} ; cp {input}.tbi {output.tbi}'
