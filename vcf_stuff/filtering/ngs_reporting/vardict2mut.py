import csv

import re
from collections import defaultdict
from os.path import join
from ngs_utils import logger
from ngs_utils.file_utils import verify_file, verify_dir, adjust_path
from ngs_utils.utils import OrderedDefaultDict

import vcf_stuff.filtering.ngs_reporting.reference_data as filt_ref_data
from vcf_stuff.filtering.ngs_reporting.utils import parse_gene_blacklists, iter_lines, parse_genes_list, check_gene_in_a_blacklist, get_anno_config


class Rule:
    def __init__(self, gene, chrom=None, start=None, end=None, length=None, ref=None,
                 required_inframe=None, indel_type=None, change=None, action=None, note=None):
        self.gene = gene
        self.chrom = chrom
        self.start = start
        self.end = end
        self.length = length
        self.ref = ref
        self.required_inframe = required_inframe
        self.indel_type = indel_type
        self.change = change
        self.action = action
        self.note = note

def parse_mut_tp53(mut_fpath):
    mut_tp53 = set()
    if verify_file(mut_fpath):
        with open(mut_fpath) as f:
            for l in f:
                l = l.strip()
                if not l:
                    continue
                line = l.split('\t')
                if not line[19] or 'p.' not in line[19]:
                    continue
                prot = line[19].replace('p.', '')
                mut_tp53.add(prot)

    return mut_tp53

def is_hotspot_nt(chrom, pos, ref, alt, hotspot_nucleotides):
    if len(ref) > len(alt) and alt != '-':
        ref = ref[1:]
        if len(alt) > 1:
            alt = alt[1:]
        else:
            alt = '-'
    elif len(alt) > len(ref) and ref != '-':
        alt = alt[1:]
        if len(ref) > 1:
            ref = ref[1:]
        else:
            ref = '-'
    key = '-'.join([chrom, str(pos), ref, alt])
    return key in hotspot_nucleotides

def is_hotspot_prot(gene, aa_chg, hotspot_proteins):
    aa_chg = aa_chg.replace('p.', '')
    if not aa_chg: return False
    key = '-'.join([gene, aa_chg])
    return key in hotspot_proteins

stop_gain_pattern = re.compile('^[A-Z]+\d+\*')
fs_pattern = re.compile('^[A-Z]+(\d+)fs')
aa_snp_chg_pattern = re.compile('^[A-Z]\d+[A-Z]$')
ins_pattern = re.compile('ins.*\*[A-Z]+$')

def parse_specific_mutations(specific_mut_fpath):
    tier_by_specific_mutations = dict()
    tier_by_type_by_region_by_gene = defaultdict(dict)
    dependent_mutations_by_gene = defaultdict(set)  # when other mutation is required

    with open(specific_mut_fpath) as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            l = l.replace('\n', '')
            if not l:
                continue
            line = l.split('\t')
            gene = line[0].upper()
            regions = re.findall(r'\d+', line[1])
            if '-' in line[1]:
                for region_num in range(int(regions[0]) + 1, int(regions[1])):
                    regions.append(str(region_num))
            if 'intron' in line[1]:
                regions = ['intron' + region for region in regions]
            for index in range(2, len(line) - 1):
                if line[index]:
                    mut = line[index]
                    tier = index - 1
                    if 'types' in mut:
                        types = mut.split(':')[1].split(',')
                        for region in regions:
                            tier_by_type_by_region_by_gene[gene][region] = dict()
                            for type in types:
                                tier_by_type_by_region_by_gene[gene][region][type] = tier
                    else:
                        mutations = []
                        if 'codon' in mut:
                            codons = re.findall(r'\d+', mut)
                            if '-' in mut and len(codons) == 2:
                                codons = range(int(codons[0]), int(codons[1]) + 1)
                            for region in regions:
                                for codon in codons:
                                    tier_by_specific_mutations['-'.join([gene, region, str(codon)])] = tier
                                    mutations.append('-'.join([gene, region, str(codon)]))
                        elif 'sens' in mut or 'res' in mut:
                            pattern = re.compile('\((\D+)\s+\D+\)')
                            sensitization = re.findall(pattern, mut)[0]  # like TKI
                            prot_chg = mut.split()[0].strip().replace('p.', '')
                            mutations = ['-'.join([gene, prot_chg])]
                            tier_by_specific_mutations['-'.join([gene, prot_chg])] = tier
                            dependent_mutations_by_gene[gene].add((sensitization, 'sens' if 'sens' in mut else 'res'))
                        else:
                            prot_chg = line[index].replace('p.', '').strip()
                            mutations = ['-'.join([gene, prot_chg])]
                            tier_by_specific_mutations['-'.join([gene, mut])] = tier

    return tier_by_specific_mutations, tier_by_type_by_region_by_gene, dependent_mutations_by_gene

class VarDict2Mut:
    def __init__(self, genome, filt_cnf, tricky_regions_dir, transcripts_fpath, reg_exp_sample=None, platform=None):
        self.all_reject_counter = OrderedDefaultDict(int)
        self.all_counter = OrderedDefaultDict(int)
        self.gene_blacklist_counter = OrderedDefaultDict(int)
        self.region_blacklist_counter = OrderedDefaultDict(int)


        compendia_fpath           = verify_file(filt_ref_data.compendia(genome), 'compendia_ms7_hotspot')
        actionable_fpath          = verify_file(filt_ref_data.actionable(genome), 'actionable')
        filter_common_snp_fpath   = verify_file(filt_ref_data.common_snp(genome), 'filter_common_snp')
        filter_common_arti_fpath  = verify_file(filt_ref_data.common_art(genome), 'filter_common_artifacts')
        splice_fpath              = verify_file(filt_ref_data.splice(genome), 'splice')
        suppressors_fpath         = verify_file(filt_ref_data.suppressors(), 'suppressors')
        oncogenes_fpath           = verify_file(filt_ref_data.oncogenes(), 'oncogenes')
        ruledir                   = verify_dir (filt_ref_data.ruledir(), 'ruledir')
        snpeffect_polymorph_fpath = verify_file(filt_ref_data.snpeffect_export_polymorphic(), 'snpeffect_export_polymorphic')
        actionable_hotspot_fpath  = verify_file(filt_ref_data.actionable_hotspot(), 'actionable_hotspot')
        specific_mutations_fpath  = verify_file(filt_ref_data.specific_mutations(), 'specific_mutations')
        last_critical_aa_fpath    = verify_file(filt_ref_data.last_critical_aa(), 'last_critical_aa')
        incidentalome_dir         = verify_dir (filt_ref_data.incidentalome_dir(), 'incidentalome')
        comments_fpath            = verify_file(filt_ref_data.ngs_reports_comments(), 'ngs_reports_comments')
        if not all([compendia_fpath          ,
                    actionable_fpath         ,
                    filter_common_snp_fpath  ,
                    filter_common_arti_fpath ,
                    splice_fpath             ,
                    suppressors_fpath        ,
                    oncogenes_fpath          ,
                    ruledir                  ,
                    snpeffect_polymorph_fpath,
                    actionable_hotspot_fpath ,
                    specific_mutations_fpath ,
                    last_critical_aa_fpath   ,
                    incidentalome_dir        ,
                    comments_fpath           ,
                    ]):
            logger.err('Error: some of the required files are not found or empty (see above)')

        self.suppressors = parse_genes_list(adjust_path(suppressors_fpath))
        self.oncogenes = parse_genes_list(adjust_path(oncogenes_fpath))

        self.reg_exp_sample = reg_exp_sample
        self.platform = platform

        transcripts_fpath = verify_file(transcripts_fpath, silent=True)
        if transcripts_fpath:
            logger.info('Using canonical transcripts from ' + transcripts_fpath)
            with open(transcripts_fpath) as f:
                self.transcripts = [tr.strip().split('.')[0] for tr in f]

        self.max_ratio = filt_cnf['max_ratio']
        self.max_sample_cnt = filt_cnf['max_sample_cnt']

        self.min_freq = filt_cnf['min_freq']  # for all variants
        self.act_min_freq = filt_cnf['act_min_freq']
        self.act_min_freq = self.act_min_freq or self.min_freq // 2
        self.germline_min_freq = filt_cnf['germline_min_freq']

        self.filt_depth = filt_cnf['filt_depth']
        self.min_vd = filt_cnf['min_vd']
        self.min_gmaf = filt_cnf['min_gmaf']

        self.keep_utr_intronic = filt_cnf['keep_utr_intronic']
        self.keep_whole_genome = filt_cnf['keep_whole_genome']
        self.keep_hla = filt_cnf['keep_hla']
        self.damage_p_value = filt_cnf.get('damage_p_value')

        logger.info('Parsing filtering data...')
        self.tp53_groups = {'Group 1': parse_mut_tp53(join(ruledir, 'DNE.txt')),
                            'Group 2': parse_mut_tp53(join(ruledir, 'TA0-25.txt')),
                            'Group 3': parse_mut_tp53(join(ruledir, 'TA25-50_SOM_10x.txt'))}

        self.splice_positions_by_gene = defaultdict(set)
        for l in iter_lines(splice_fpath):
            pos, g = l.split('\t')
            self.splice_positions_by_gene[g].add(pos)

        self.last_critical_aa_pos_by_gene = dict()
        for l in iter_lines(last_critical_aa_fpath):
            g, aa_pos, _ = l.split('\t')
            self.last_critical_aa_pos_by_gene[g] = int(aa_pos)

        self.filter_snp = set()
        for l in iter_lines(filter_common_snp_fpath):
            fields = l.split('\t')
            self.filter_snp.add('-'.join(fields[1:5]))

        self.snpeff_snp = set()
        self.snpeff_snp_rsids = set()
        for l in iter_lines(snpeffect_polymorph_fpath):
            fields = l.split('\t')
            snpeff_aachg = fields[2]
            snpeff_rsid = fields[5]
            if len(fields) > 11 and fields[11]:
                snpeff_gene = fields[11]
                self.snpeff_snp.add('-'.join([snpeff_gene, snpeff_aachg]))
            elif snpeff_rsid != '-':
                self.snpeff_snp_rsids.add(snpeff_rsid)

        self.filter_artifacts = set()
        self.filter_rules_by_gene = defaultdict(list)
        for l in iter_lines(filter_common_arti_fpath):
            fields = l.split('\t')
            if fields[5] == 'rule':
                gene, chrom, start, end, action, _, _, _, note  = fields[:9]
                rule = Rule(gene, chrom=chrom, start=int(start), end=int(end), action=action, note=note)
                self.filter_rules_by_gene[gene].append(rule)
            else:
                gene, chrom, start, ref, alt = fields[:5]
                self.filter_artifacts.add('-'.join([chrom, start, ref, alt]))

        self.actionable_hotspot_by_gene = defaultdict(dict)
        self.common_snps_by_gene = defaultdict(set)
        with open(actionable_hotspot_fpath) as f:
            for l in f:
                l = l.replace('\n', '')
                if not l or l.startswith('##'):
                    continue
                fields = l.split('\t')
                gene = fields[0]
                prot_change = fields[1]
                if gene.startswith('#'):  # VUS, No special treatment for now
                    gene = gene[1:]
                elif gene.startswith('^'):
                    gene = gene[1:]
                    self.common_snps_by_gene[gene].add(prot_change)
                else:
                    is_somatic = fields[2] == 'somatic'
                    self.actionable_hotspot_by_gene[gene][prot_change] = 'somatic' if is_somatic else 'germline'

        self.ngs_reports_comments = defaultdict(dict)
        with open(comments_fpath) as f:
            for r in csv.DictReader((row for row in f if not row.startswith('#')), delimiter='\t'):
                gene = r['Gene']
                prot_change = r['AA_Change']
                if gene.startswith('^'):
                    gene = gene[1:]  # remove leading ^ character, e.g. ^EGFR -> EGFR
                    is_somatic = 'somatic' in r['Note']
                    self.actionable_hotspot_by_gene[gene][prot_change] = 'somatic' if is_somatic else 'germline'
                else:
                    self.ngs_reports_comments[gene][prot_change] = r['Note']

        self.act_somatic = dict()
        self.act_germline = set()
        self.rules = defaultdict(list)
        for l in iter_lines(actionable_fpath):
            fields = l.split('\t')

            if fields[7] == 'germline':
                key = '-'.join(fields[1:5])
                self.act_germline.add(key)

            elif fields[7] == 'somatic':
                change = fields[8].strip()
                if fields[6] == 'rule':
                    if fields[4] == '*' and len(fields[3]) == 1:
                        key = '-'.join(fields[1:4])
                        self.act_somatic[key] = change
                    else:
                        indel_type = ''
                        if 'indel' in fields[5]: indel_type = 'indel'
                        elif 'ins' in fields[5]: indel_type = 'ins'
                        elif 'del' in fields[5]: indel_type = 'del'
                        rule = Rule(gene=fields[0],
                                    chrom=fields[1],
                                    start=int(fields[2]),
                                    end=int(fields[3]),
                                    length=int(fields[4]),
                                    required_inframe='inframe' in fields[5],
                                    indel_type=indel_type,
                                    change=change)
                        self.rules[rule.gene].append(rule)
                    # elif fields[5] == inframe_del:
                    #     self.rules[inframe_del].setdefault(fields[0], []).append([fields[1]] + [int (f) for f in fields[2:5]])
                    # elif fields[5] == inframe_ins:
                    #     self.rules[inframe_ins].setdefault(fields[0], []).append([fields[1]] + [int (f) for f in fields[2:5]])

                else:
                    key = '-'.join(fields[1:5])
                    self.act_somatic[key] = change

        self.hotspot_nucleotides = set()
        self.hotspot_proteins = set()
        for l in iter_lines(compendia_fpath):
            fields = l.split('\t')
            if fields[5].startswith('g.'):
                continue
            self.hotspot_nucleotides.add('-'.join(fields[1:5]))
            if not fields[6]:
                continue
            self.hotspot_proteins.add('-'.join([fields[0], fields[6]]))

        logger.info('Parsing gene blacklists...')
        anno_cfg = get_anno_config()
        self.gene_blacklists_by_reason = parse_gene_blacklists(anno_cfg['blacklist']['genes'], incidentalome_dir)
        for r in self.gene_blacklists_by_reason.keys():
            self.gene_blacklist_counter[r] = 0
        self.gene_blacklist_counter['hardfilter'] = 0
        # self.gene_to_soft_filter = list(iter_lines(join(incidentalome_dir, 'soft_filter.txt')))

        # self.region_blacklists_by_reason = dict()
        # if tricky_regions_dir:
        #     info('Parsing region blacklists...')
        #     self.region_blacklists_by_reason = load_tricky_regions(anno_cfg['blacklist']['regions'], tricky_regions_dir)
        #     for r in self.region_blacklists_by_reason.keys():
        #         self.region_blacklist_counter[r] = 0

        logger.info('Parsing actionable rules and specific mutations...')
        self.tier_by_specific_mutations, self.tier_by_type_by_region_by_gene, self.sensitizations_by_gene\
            = parse_specific_mutations(specific_mutations_fpath)

        if not all([self.rules, self.splice_positions_by_gene, self.act_somatic, self.act_germline, self.actionable_hotspot_by_gene]):
            if not self.rules:
                logger.err('No rules, cannot proceed')
            if not self.splice_positions_by_gene:
                logger.err('No tp53_positions, cannot proceed')
            if not self.act_somatic:
                logger.err('No act_somatic, cannot proceed')
            if not self.act_germline:
                logger.err('No act_germline, cannot proceed')
            if not self.actionable_hotspot_by_gene:
                logger.err('No actionable_hotspots, cannot proceed')

        self.status = None
        self.reason_by_status = None

        self.output_f = None
        self.fm_output_f = None
        self.rejected_output_f = None

    aa_chg_trim_pattern = re.compile('^([A-Z]\d+)[A-Z?]$')
    def check_actionable(self, chrom, pos, ref, alt, gene, aa_chg, cosm_aa_chg, af, clnsig):
        change_len = len(alt) - len(ref)

        key = '-'.join([chrom, str(pos), ref, alt])
        if key in self.act_somatic:
            return 'act_somatic'
        if key in self.act_germline and af >= self.germline_min_freq:
            return 'act_germline'

        if len(ref) == 1 and change_len == 0:  # SNP
            key = '-'.join([chrom, str(pos), ref])
            if key in self.act_somatic:
                return 'act_somatic'

        if gene in self.actionable_hotspot_by_gene and \
                (VarDict2Mut.aa_chg_trim_pattern.match(aa_chg) or re.compile('^(M1)\?$').match(aa_chg)):
            act_hotspot_by_aa_chg = self.actionable_hotspot_by_gene[gene]
            status = act_hotspot_by_aa_chg.get(aa_chg)
            if status is not None:
                if status == 'somatic' and af > self.act_min_freq:
                    return 'act_hotspot_somatic'
                elif status == 'germline' and af > self.germline_min_freq:
                    return 'act_hotspot_germline'
            aa_chg_trim = re.findall(VarDict2Mut.aa_chg_trim_pattern, aa_chg)[0]
            status = act_hotspot_by_aa_chg.get(aa_chg_trim)
            if status is not None:
                return 'act_hotspot_' + status

        if gene == 'TP53':
            tp53_group = self.classify_tp53(aa_chg, pos, ref, alt)
            if tp53_group is not None:
                return 'act_somatic_tp53_group_' + str(tp53_group)

        if gene in self.rules:
            for r in self.rules[gene]:
                if change_len >= r.length and r.start <= pos + len(ref) - 1 and pos <= r.end:
                    if r.required_inframe and change_len % 3 != 0:
                        continue
                    if any([r.indel_type == 'ins' and change_len > 0,
                            r.indel_type == 'del' and change_len < 0,
                            r.indel_type == 'indel' and change_len != 0]):
                        return 'somatic'
        return None

    def classify_tp53(self, aa_chg, pos, ref, alt):
        aa_chg = aa_chg.replace(' ', '')
        if str(pos) in self.splice_positions_by_gene['TP53'] and len(ref) == 1 and len(alt) == 1:
            return 6
        aa_chg = aa_chg.replace('p.', '')
        aa_num = 0
        if aa_chg:
            aa_num_str = re.sub('[^0-9]', '', aa_chg)
            if not aa_num_str:
                logger.err('TP53: cannot parse aa num from aa_chg=' + str(aa_chg))
            else:
                aa_num = int(aa_num_str)
        if aa_snp_chg_pattern.match(aa_chg):
            for i in [1, 2, 3]:
                if aa_chg in self.tp53_groups['Group ' + str(i)]:
                    return i
        elif stop_gain_pattern.match(aa_chg):
            if aa_num < 359:
                return 4
        elif fs_pattern.match(aa_chg):
            if aa_num < 359:
                return 5
        return None

    def check_rob_hedley_actionable(self, gene, aa_chg, effect, region):
        if aa_chg:
            gene_aachg = '-'.join([gene, aa_chg])
            if gene_aachg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_aachg]
                return 'actionable' if tier == 1 else 'tier2'

        if region and effect in ['HIGH', 'MODERATE']:
            codon = re.sub('[^0-9]', '', aa_chg)
            gene_codon_chg = '-'.join([gene, region, codon])
            if gene_codon_chg in self.tier_by_specific_mutations:
                tier = self.tier_by_specific_mutations[gene_codon_chg]
                return ('act' if tier == 1 else 'tier2') + '_codon_' + codon + '_in_exon_' + region

    def check_by_type_and_region(self, cdna_chg, region, gene):
        types_by_region = self.tier_by_type_by_region_by_gene.get(gene)
        if types_by_region:
            for type_ in types_by_region.get(region, []):
                if type_ in cdna_chg:
                    tier = types_by_region[region][type_]
                    return ('act' if tier == 1 else 'tier2') + '_' + type_ + '_in_gene_' + gene
        return False

    def _check_artifacts(self, chrom, pos, id_field, gene, aa_chg):
        snps = re.findall(r'rs\d+', id_field)
        if any(snp in self.snpeff_snp_rsids for snp in snps):
            return 'snp in snpeffect_export_polymorphic'

        for r in self.filter_rules_by_gene.get(gene, []):
            if r.action == 'ignore' and chrom == r.chrom and r.start <= pos <= r.end:
                if r.note:
                    return 'filter artifacts: ' + r.note
                else:
                    return 'filter artifacts: ' + r.action

        if gene in self.common_snps_by_gene and aa_chg in self.common_snps_by_gene[gene]:
            return 'common SNP'

        if gene in self.ngs_reports_comments and aa_chg in self.ngs_reports_comments[gene]:
            return 'ngs report comment: ' + self.ngs_reports_comments[gene][aa_chg]

        return None

    def print_mutation(self, status, reasons, blacklisted_reasons, fields, fm_data):
        self.apply_counter('lines_written')
        self.apply_counter(status)

        if fm_data and self.fm_output_f:
            sample, platform, gene, pos, cosm_aa_chg, aa_chg, cdna_chg, chrom, depth, allele_freq = fm_data
            self.fm_output_f.write('\t'.join([sample, platform, 'short-variant', gene, status, aa_chg, cdna_chg,
                                   chrom + ':' + pos, str(depth), str(allele_freq * 100),
                                   '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-'])  + '\n')
        self.output_f.write('\t'.join(fields + [status] + [', '.join(reasons)] + [', '.join(blacklisted_reasons)]) + '\n')

    def apply_counter(self, reason):
        self.all_counter[reason] += 1

    def reject_mutation(self, reason, fields):
        self.apply_reject_counter(reason)
        if self.rejected_output_f:
            self.rejected_output_f.write('\t'.join(fields + ['rejected'] + [reason] + ['']) + '\n')

    def apply_reject_counter(self, reason):
        self.all_reject_counter[reason] += 1

    def apply_gene_blacklist_counter(self, reason):
        self.gene_blacklist_counter[reason] += 1

    def apply_region_blacklist_counter(self, reason):
        self.region_blacklist_counter[reason] += 1

    def check_blacklist_genes(self, gene_name, aa_pos=None):
        reasons = []
        for reason, data in self.gene_blacklists_by_reason.items():
            if check_gene_in_a_blacklist(gene_name, data, aa_pos):
                reasons.append(reason)
        return reasons

    def filter(self, rec, tumor_indices):
        id_field = rec.ID or ""
        chrom = rec.CHROM
        pos = rec.POS
        gene = rec.INFO.get('PCGR_SYMBOL')
        aa_chg = ''
        ref = rec.REF
        alt = rec.ALT[0]
        af = rec.INFO['TUMOR_AF']
        clnsig = rec.INFO.get('PCGR_CLINVAR_CLNSIG')

        if 'chr' not in chrom: chrom = 'chr' + chrom
        nt_chg_key = '-'.join([chrom, str(pos), ref, alt])


        fail_reason = self._check_artifacts(chrom, pos, id_field, gene, aa_chg)
        if fail_reason:
            rec.INFO['AZ_artefact'] = fail_reason

        # cosmic_counts = map(int, fields[cosmcnt_col].split()) if cosmcnt_col is not None else None
        # cosm_aa_chg = self.check_by_var_class(var_class, cosm_aa_chg, cosmic_counts)

        # aa_chg = self.check_by_effect(var_type, aa_chg, cdna_chg, effect)
        # aa_pos = None
        # if aa_chg_pos_regexp.match(aa_chg):
        #     aa_pos = int(aa_chg_pos_regexp.findall(aa_chg)[0])

        if is_hotspot_nt(chrom, pos, ref, alt, self.hotspot_nucleotides):
            print('Setting AZ_hotspot')
            rec.INFO['AZ_hotspot'] = 'hotspot_nucl_change'
            print('Set AZ_hotspot')
        # elif is_hotspot_prot(gene, aa_chg, self.hotspot_proteins):
        #     self.update_status('likely', 'hotspot_AA_change')

        #################################
        # Checking actionable mutations #
        #################################
        actionable_status = \
            self.check_actionable(chrom, pos, ref, alt, gene, cosm_aa_chg=None, aa_chg=aa_chg, af=af, clnsig=clnsig) or \
            self.check_rob_hedley_actionable(gene, aa_chg, '', '') or \
            self.check_by_type_and_region('', '', gene)

        if actionable_status:
            if 'germline' in actionable_status and af < self.germline_min_freq:
                rec.INFO['AZ_artefact'] = 'act germline and AF < ' + str(self.germline_min_freq)
            elif af < self.act_min_freq:
                rec.INFO['AZ_artefact'] = 'act somatic and AF < ' + str(self.act_min_freq)
            else:
                rec.INFO['AZ_hotspot'] = actionable_status

        else:
            if nt_chg_key in self.filter_snp:
                rec.INFO['AZ_artefact'] = 'not act and in filter_common_snp'
            if nt_chg_key in self.filter_artifacts and af < 0.35:
                rec.INFO['AZ_artefact'] = 'not act and in filter_artifacts and AF < 0.35'
            # gmaf = fields[headers.index('GMAF')]
            # if gmaf and all(not g or float(g) > self.min_gmaf for g in gmaf.split(',')):
            #     self.reject_mutation('not act and all GMAF > ' + str(self.min_gmaf), fields)
            #     continue
            # if var_class == 'dbSNP':
            #     self.reject_mutation('clnsig dbSNP', fields)
            #     continue
            # if '-'.join([gene, aa_chg]) in self.snpeff_snp and var_class != 'ClnSNP_known':
            #     self.reject_mutation('not known and in SnpEff SNP database', fields)
            #     continue

            # if self.min_freq and af < self.min_freq:
            #     self.reject_mutation('not act and AF < ' + str(self.min_freq) + ' (min_freq)', fields)
            #     continue

            if self.check_blacklist_genes(gene):
                rec.INFO['AZ_artefact'] = 'not known and blacklist gene'

        # Ignore any variants that occur after last known critical amino acid
        # if aa_pos is not None:
        #     if gene in self.last_critical_aa_pos_by_gene and aa_pos >= self.last_critical_aa_pos_by_gene[gene]:
        #         self.reject_mutation('variants occurs after last known critical amino acid', fields)
        #         continue

        # blacklisted_reasons = self.check_callability_regions(chrom=chrom, start=pos - 1, end=pos - 1 + len(ref)

        return rec
