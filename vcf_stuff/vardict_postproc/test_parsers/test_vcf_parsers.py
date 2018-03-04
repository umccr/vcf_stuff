import sys


def proc_fields(ref, alt, af, msi):
    # Filter low AF MSI
    if msi:
        msi = float(msi)
        change_len = abs(len(ref) - len(alt))
        if change_len == 1 and msi > 1:
            msi_fail = any([
                msi <=  2 and af < 0.005,
                msi <=  4 and af < 0.01,
                msi <=  7 and af < 0.03,
                msi ==  8 and af < 0.06,
                msi ==  9 and af < 0.125,
                msi == 10 and af < 0.175,
                msi == 11 and af < 0.25,
                msi == 12 and af < 0.3,
                msi >  12 and af < 0.35])
            if msi_fail:
                return False
        elif change_len == 3 and msi >= 5 and af < 0.1:  # ignore low AF in 3nt MSI region
            return False
    return True


def use_pyvcf(vcf_file):
    """ Working, but doesn't write header. Very slow (first test was >17m).
    """
    import vcf
    import gzip
    f = gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file)
    with f as f:
        vcf_reader = vcf.Reader(f)
        vcf_writer = vcf.Writer(sys.stdout, vcf_reader)
        for rec in vcf_reader:
            msi_fail = proc_fields(rec.REF, rec.ALT[0], rec.samples[0]['AF'], rec.INFO['MSI'])
            if msi_fail:
                rec.FILTER.append('MSI_FAIL')
            vcf_writer.write_record(rec)

def use_cyvcf(vcf_file):
    """ Not installing (py3 not supported, py2 installs but import doesn't work)
    """
    import cyvcf  # need to reinstall instead of pyvcf
    import gzip
    f = gzip.open(vcf_file) if vcf_file.endswith('.gz') else open(vcf_file)
    with f as f:
        vcf_reader = cyvcf.Reader(f)
        vcf_writer = cyvcf.Writer(sys.stdout, vcf_reader)
        for rec in vcf_reader:
            msi_fail = proc_fields(rec.REF, rec.ALT[0], rec.samples[0]['AF'], rec.INFO['MSI'])
            if msi_fail:
                rec.FILTER.append('MSI_FAIL')
            vcf_writer.write_record(rec)

def use_pysam(vcf_file):
    """ Working.
        Time:           3:34.04
    """
    import pysam
    vcf = pysam.VariantFile(vcf_file)
    vcf.header.filters.add('MSI_FAIL', None, None, '')
    sys.stdout.write(str(vcf.header))
    for rec in vcf:
        msi_fail = proc_fields(rec.ref, rec.alts[0], rec.samples.values()[0]['AF'], rec.info['MSI'])
        if msi_fail:
            rec.filter.add('MSI_FAIL')
        sys.stdout.write(str(rec))

def use_cyvcf2(vcf_file, vcf_out=None):
    """ Working.
        File out:       2:17.51
        stdout + bgzip: 2:50.35
    """
    from cyvcf2 import VCF, Writer

    vcf = VCF(vcf_file)
    vcf.add_filter_to_header({'ID': 'MSI_FAIL', 'Description': 'Possible homopolymer artefact'})
    if vcf_out:
        w = Writer(vcf_out, vcf)
    else:
        w = None
        sys.stdout.write(vcf.raw_header)
    for rec in vcf:
        msi_fail = proc_fields(rec.REF, rec.ALT[0], rec.format('AF')[0][0], rec.INFO['MSI'])
        if msi_fail:
            filters = rec.FILTER.split(';') if rec.FILTER else []
            filters.append('MSI_FAIL')
            rec.FILTER = ';'.join(filters)
        if w:
            w.write_record(rec)
        else:
            sys.stdout.write(str(rec))
    if w:
        w.close()

def proc_line(line):
    import re
    if line.startswith('#'):
        return line
    else:
        c, p, i, ref, alt, q, filt, info, frmt = line.split('\t')[:9]
        samples = line.split('\t')[9:]
        m = re.match(r'.*;?MSI=(?P<msi>\d)?.*', info)
        if m:
            msi = m.group('msi')
        else:
            return line
        af_idx = frmt.split(':').index('AF')
        try:
            af = float(samples[0].split(':')[af_idx])
        except ValueError:
            sys.stderr.write('Could not parse AF from line ' + line + '\n')
            sys.exit(1)
        msi_fail = proc_fields(ref, alt, af, msi)
        if msi_fail:
            filt = ';'.join(filt.split(';') + [';MSI_FAIL'])
        return '\t'.join([c, p, i, ref, alt, q, filt, info, frmt] + samples)

def use_python(vcf_file):
    """ Working.
        Time:              2:26.87
        Through pythonpy:  2:33.13
    """
    import gzip
    f = gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file)
    with f as f:
        for line in f:
            if line.startswith('#CHROM'):
                sys.stdout.write('##FILTER=<ID=MSI_FAIL,Description="Possible homopolymer artefact">\n')
            line = proc_line(line)
            sys.stdout.write(line)


def main():
    vcf_file = sys.argv[1]
    cmd = sys.argv[2]

    if cmd == 'cyvcf':
        use_cyvcf(vcf_file)
    if cmd == 'python':
        use_python(vcf_file)
    if cmd == 'pysam':
        use_pysam(vcf_file)
    if cmd == 'pyvcf':
        use_pyvcf(vcf_file)
    if cmd == 'cyvcf2':
        vcf_out = sys.argv[3] if len(sys.argv) > 3 else None
        use_cyvcf2(vcf_file, vcf_out)


if __name__ == '__main__':
    main()
