#!/usr/bin/env python

import sys
from pprint import pprint
import click
from os.path import isfile, join
from cyvcf2 import VCF, Writer
import re
import numpy as np
from umccrise.utils import get_sample_ids


@click.command()
@click.argument('input_file', type=click.Path(exists=True))
@click.option('-o', 'output_file', type=click.Path())
def main(input_file, output_file=None):
    vcf = VCF(input_file, gts012=True)

    if output_file:
        w = Writer(output_file, vcf)
    else:
        w = None

    tumor_index, control_index = get_sample_ids(input_file)

    tumor_prefix = 'TUMOR_'
    normal_prefix = 'NORMAL_'

    # Add headers
    new_header = []
    for h in vcf.raw_header.split('\n'):
        if h.startswith('#CHROM'):
            for tag in ['AF', 'DP', 'MQ']:
                type = 'Integer' if tag == 'DP' else 'Float'
                tumor_tag = tumor_prefix + tag
                normal_tag = normal_prefix + tag

                new_header.append(f'##INFO=<ID={tumor_tag},Number=1,Type={type},Description="{tag} in tumor sample">')
                # Add headers to original cyvcf2 header just to avoid errors
                vcf.add_info_to_header({'ID': tumor_tag, 'Description': '', 'Type': type, 'Number': '1'})
                if control_index is not None:
                    new_header.append(f'##INFO=<ID={normal_tag},Number=1,Type={type},Description="{tag} in control sample">')
                    # Add headers to original cyvcf2 header just to avoid errors
                    vcf.add_info_to_header({'ID': normal_tag, 'Description': '', 'Type': type, 'Number': '1'})

        if not any(h.startswith(f'##INFO=<ID={p+tag},') for tag in ['AF', 'DP', 'MQ'] for p in [tumor_prefix, normal_prefix]):
            new_header.append(h)

    sys.stdout.write('\n'.join(new_header))

    # Go through each record and add new INFO fields
    for rec in vcf:
        af, dp, mq = _collect_vals_per_sample(rec, control_index, tumor_index)

        for t, v in zip(['AF', 'DP', 'MQ'], [af, dp, mq]):
            rec.INFO[tumor_prefix + t] = str(v[tumor_index])
            if control_index is not None:
                if len(v) <= control_index:
                    sys.stderr.write(f'Warning: for tag {t}, len of v={len(v)} is less than index {control_index} of control sample. Record {v}\n')
                elif v[control_index] is not None:
                    rec.INFO[normal_prefix + t] = str(v[control_index])

        if w:
            w.write_record(rec)
        else:
            sys.stdout.write(str(rec))

    if w:
        w.close()
    # vcf.close()


''' 
VarDict
FORMAT/DP,       FORMAT/AF,                                FORMAT/MQ (somatic), INFO/MQ (germline)
                     
Mutect2                     
FORMAT/DP,       FORMAT/AF,                                FORMAT/MMQ
                     
Freebayes                                               
FORMAT/DP        FORMAT/AD = ref_count,alt_count           INFO/MQM+INFO/MQMR
                     
GATK-Haplotype                                          
FORMAT/DP        FORMAT/AD = ref_count,alt_count           FORMAT/MMQ
         
Strelka2 - germline                      
SNV:     
sum(alt_counts)  FORMAT/AD = ref_count,alt_counts          INFO/MQ
INDEL:     
sum(alt_counts)  FORMAT/AD = ref_count,alt_counts          INFO/MQ

Strelka2 - somatic                 
SNV:
FORMAT/DP        FORMAT/{ALT}U[0] = alt_count(tier1,tier2) INFO/MQ
INDEL:
FORMAT/DP        FORMAT/TIR = alt_count(tier1,tier2)       INFO/MQ
'''


def _collect_vals_per_sample(rec, control_index, tumor_index):
    dp = None
    af = None
    mq = None

    # If FORMAT/AF exists, report it as af. Else, check FORMAT/AD. If not, check FORMAT/*U
    if 'AF' in rec.FORMAT:
        af = rec.format('AF')[:,0]
        dp = rec.format('DP')[:,0]
    else:
        # strelka2 germline?
        if 'AD' in rec.FORMAT:
            alt_counts = rec.format('AD')[:,1]  # AD=REF,ALT so 1 is the position of ALT
            dp = np.sum(rec.format('AD')[:,0:], axis=1)
        # strelka2 somatic?
        else:
            if rec.is_snp:
                alt_counts = rec.format(rec.ALT[0] + 'U')[:,0]
            else:
                alt_counts = rec.format('TIR')[:,0]
            dp = rec.format('DP')[:,0]
        af = np.true_divide(alt_counts, dp, out=np.zeros(alt_counts.shape), where=dp!=0)

    try:
        mq = rec.format('MQ', float)[:,0]  # VarDict has an incorrect MQ header (with Integer type instead of Float), so need to specify "float" type here explicitly otherwise MQ won't be parsed
    except:
        try:
            mq = rec.format('MMQ')[:,0]
        except:
            mq = [None for _ in [control_index, tumor_index]]
            mq[tumor_index] = rec.INFO.get('MQ')

    return af, dp, mq


def _parse_tag(rec, header_by_tag, tag, d, tumor_index, control_index):
    header = header_by_tag[tag]
    if header:
        if header['HeaderType'] == 'FORMAT':
            data = rec.format(header['ID'], header['python_type'])
            d[tag]['tumor'] = str(data[tumor_index][0])

            sample_dim = data.shape[0]
            if sample_dim >= 2:
                d[tag]['normal'] = str(data[control_index][0])

            if d[tag]['tumor']  == '-2147483648': d[tag]['tumor']  = -1
            if d[tag]['normal'] == '-2147483648': d[tag]['normal'] = -1
        else:
            d[tag]['tumor'] = str(rec.INFO[header['ID']])


if __name__ == '__main__':
    main()
