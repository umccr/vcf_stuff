import subprocess
import sys
from cyvcf2 import VCF, Writer
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import get_ungz_gz


def count_vars(vcf_path, filter=None):
    cmd = f'bcftools view -H {f"-f {filter} " if filter else " "}{vcf_path} | wc -l'
    return int(subprocess.check_output(cmd, shell=True).strip())


def vcf_contains_field(vcf_path, field, col=None):
    # col is FILTER, FORMAT, INFO, or None (=any of three)
    if col is None and '/' in field:
        col, field = field.split('/')
    if col is not None:
        return f'##{col}=<ID={field},' in VCF(vcf_path).raw_header
    return f'=<ID={field},' in VCF(vcf_path).raw_header


def iter_vcf(input_file, output_file, proc_rec, proc_hdr=None, postproc_hdr=None, **kwargs):
    """
    :param input_file: path to input VCF file
    :param output_file: path to output VCF file (can be .vcf or .vcf.gz, but it will always bgzip/tabix and write with .vcf.gz extention)
    :param proc_rec: a function to process a single cyvcf Record object. Returns either a (new) Record object to write, or None to indicate that the record should be discarded
    :param proc_hdr: a function to process cyvcf object once (i.e. to add values to the header with vcf.add_info_to_header, etc)
    :param postproc_hdr: a function to postprocess finalized header string (vcf.rawheader), e.g. in order to remove values
    :param kwargs: any paramters to pass directly into proc_rec
    """
    vcf = VCF(input_file, gts012=True)
    if proc_hdr is not None:
        proc_hdr(vcf)

    # w = None
    if output_file is not None:
        out_ungz, out_gz = get_ungz_gz(output_file)
        # w = Writer(out_ungz, vcf)
        # w.write_header()
        w = open(out_ungz, 'w')
    else:
        # sys.stdout.write(vcf.raw_header)
        w = sys.stdout

    header = vcf.raw_header
    if postproc_hdr is not None:
        header = postproc_hdr(header)
    w.write(header)

    for rec in vcf:
        if proc_rec:
            rec_res = proc_rec(rec, vcf, **kwargs)
            if rec_res is not None:
                # if w is not None:
                #     sys.stderr.write('Writing record', rec_res, '\n')
                #     w.write_record(rec_res)
                # else:
                #     print(rec_res)
                # sys.stderr.write(f'Writing record {rec_res}\n')
                w.write(f'{rec_res}')

    sys.stderr.write(f'Finished writing {output_file}\n')
    vcf.close()
    if output_file is not None:
        w.close()
        run_simple(f'bgzip -f {out_ungz} && tabix -f -p vcf {out_gz}')
        sys.stderr.write(f'Compressed {output_file}\n')


def iter_vcf__pysam(input_file, proc_rec=None, proc_hdr=None, output_file=None):
    import pysam
    import sys

    vcf = pysam.VariantFile(input_file)
    if output_file:
        w = open(output_file, 'w')
    else:
        w = sys.stdout

    # Header
    if proc_hdr is not None:
        proc_hdr(vcf)
    w.write(str(vcf.header))

    # Records
    for rec in vcf:
        if proc_rec:
            rec_res = proc_rec(rec)
            if rec_res is not None:
                print(rec_res)
                w.write(str(rec_res))

    vcf.close()

    if output_file:
        w.close()
        out_ungz, out_gz = get_ungz_gz(output_file)
        run_simple(f'bgzip -f {out_ungz} && tabix -f -p vcf {out_gz}')

