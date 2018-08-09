
def make_normalise_cmd(input_file, output_file, reference_fasta):
    return (
        f'bcftools norm -m \'-\' {input_file} -Ov -f {reference_fasta}'     # split multiallelic ALT and left-aligns indels
        f' | vcfallelicprimitives -t DECOMPOSED --keep-geno --keep-info'    # split MNP into single SNPs
        f' | vcfstreamsort'
        f' | grep -v "##INFO=<ID=TYPE,Number=1"'
        f' | bgzip -c > {output_file}'
        f' && tabix -p vcf {output_file}')


