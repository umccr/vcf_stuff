from os.path import isfile, join, basename, splitext
import shutil
from ngs_utils import logger as log
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.vcf_utils import get_sample_ids, get_sample_names, iter_vcf
from reference_data import api as refdata
import cyvcf2
from vcf_stuff.filtering import add_cyvcf2_filter, package_path


localrules: sage


SAMPLE = config['sample']
GENOME = config['genome']
TUMOR_BAM = config['tumor_bam']
NORMAL_BAM = config['normal_bam']
EXISTING_VCF = config.get('existing_vcf', None)
OUTPUT_EXISTING_SAGED_VCF = config['output_existing_saged_vcf']
OUTPUT_SAGE_VCF = config['output_sage_vcf']
assert EXISTING_VCF.endswith('.vcf.gz'), EXISTING_VCF
assert OUTPUT_EXISTING_SAGED_VCF.endswith('.vcf.gz'), OUTPUT_EXISTING_SAGED_VCF
assert OUTPUT_SAGE_VCF.endswith('.vcf.gz'), OUTPUT_SAGE_VCF
if EXISTING_VCF:
    T_NAME, N_NAME = get_sample_names(
        EXISTING_VCF,
        provided_tumor_name=config.get('tumor_vcf_sample'),
        provided_normal_name=config.get('normal_vcf_sample')
    )
else:
    T_NAME = config.get('tumor_vcf_sample') or splitext(basename(TUMOR_BAM))
    N_NAME = config.get('normal_vcf_sample') or splitext(basename(NORMAL_BAM))

refdata.find_genomes_dir(config.get('input_genomes_url'))
HOTSPOTS_VCF = config.get('hotspots_vcf', refdata.get_ref_file(GENOME, key='hotspots'))
if config.get('call_inframe') is True:
    CODING_REGIONS = config.get('coding_regions', refdata.get_ref_file(GENOME, key='coding_regions'))
    log.warn(f'Calling inframe indels in {CODING_REGIONS}')
else:
    CODING_REGIONS = '<(echo "")'
    log.warn(f'Not calling inframe indels')

rule all:
    input:
        vcf = OUTPUT_EXISTING_SAGED_VCF if EXISTING_VCF else [],
        tbi = (OUTPUT_EXISTING_SAGED_VCF + '.tbi') if EXISTING_VCF else [],
        sage_vcf = OUTPUT_SAGE_VCF,
        sage_tbi = OUTPUT_SAGE_VCF + '.tbi',


rule run_sage:
    input:
        tumor_bam    = TUMOR_BAM,
        tumor_bai    = TUMOR_BAM + '.bai',
        normal_bam   = NORMAL_BAM,
        normal_bai   = NORMAL_BAM + '.bai',
        ref_fa       = refdata.get_ref_file(GENOME, key='fa'),
        hotspots_vcf = HOTSPOTS_VCF,
    output:
        sage_vcf = f'work/call/{SAMPLE}-sage.vcf.gz',
        sage_tbi = f'work/call/{SAMPLE}-sage.vcf.gz.tbi',
    params:
        jar = join(package_path(), 'sage-1.0.jar'),
        tumor_sname  = T_NAME,
        normal_sname = N_NAME,
        xms = 2000,
        xmx = 19000,
        coding_bed = CODING_REGIONS,
    resources:
        mem_mb = 20000
    group: "sage"
    shell:
        'java -Xms{params.xms}m -Xmx{params.xmx}m -jar {params.jar} '
        '-tumor {params.tumor_sname} -tumor_bam {input.tumor_bam} '
        '-reference {params.normal_sname} -reference_bam {input.normal_bam} '
        '-known_hotspots <(bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" {input.hotspots_vcf}) '
        '-coding_regions {params.coding_bed} '
        '-ref_genome {input.ref_fa} '
        '-out {output.sage_vcf} '
        '&& tabix -f -p vcf {output.sage_vcf}'


if EXISTING_VCF:
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
            idx = get_sample_ids(input.vcf, provided_t_name=T_NAME, provided_n_name=N_NAME)
            tumor_index, normal_index = idx[:2]
            assert sorted([tumor_index, normal_index]) == [0, 1]
            sample_in_order = [None, None]
            sample_in_order[tumor_index] = T_NAME
            sample_in_order[normal_index] = N_NAME
            shell(f'bcftools view -s {",".join(sample_in_order)} {input.sage_vcf} -Oz -o {output.sage_vcf} '
                  f'&& tabix -p vcf {output.sage_vcf}')


    rule sage_pass:
        input:
            sage_vcf = rules.sage_reorder_samples.output.sage_vcf,
            sage_tbi = rules.sage_reorder_samples.output.sage_tbi,
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

            def proc_rec(rec, vcf):
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
            iter_vcf(input.vcf, output.vcf, proc_rec, proc_hdr=proc_hdr)

    rule copy_result:
        input:
            vcf = rules.annotate_from_sage.output.vcf,
            tbi = rules.annotate_from_sage.output.tbi,
        output:
            vcf = OUTPUT_EXISTING_SAGED_VCF,
            tbi = OUTPUT_EXISTING_SAGED_VCF + '.tbi'
        shell:
            'cp {input.vcf} {output.vcf} ; '
            'cp {input.tbi} {output.tbi}'


rule sage:
    input:
        sage_vcf = rules.run_sage.output.sage_vcf,
        sage_tbi = rules.run_sage.output.sage_tbi,
    output:
        vcf = OUTPUT_SAGE_VCF,
        tbi = OUTPUT_SAGE_VCF + '.tbi'
    shell:
        'cp {input.sage_vcf} {output.vcf} ; '
        'cp {input.sage_tbi} {output.tbi} ;'


onsuccess:
    print("sage workflow finished! Deleting .snakemake/metadata")
    shutil.rmtree(".snakemake/metadata")

