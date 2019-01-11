import sys

from cyvcf2 import VCF, Writer

from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.vcf_utils import get_sample_ids


def iter_vcf(input_file, output_file, proc_rec, proc_hdr=None, **kwargs):
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
    w.write(vcf.raw_header)

    for rec in vcf:
        if proc_rec:
            rec_res = proc_rec(rec, **kwargs)
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

