from collections import OrderedDict

text = '''
bwa-strelka2             3206949  4038   4171    504322  30400  9861   
bwa-vardict              3195487  27840  15633   473815  16530  40368  
bwa-gatk-haplotype       3209073  8345   2047    506816  7771   7367   
bwa-ensemble             3209551  4699   1569    503714  4515   10469  
minimap2-strelka2        3202865  3342   8255    503762  28823  10421  
minimap2-vardict         3194962  26981  16158   475068  16338  39115  
minimap2-gatk-haplotype  3207496  7050   3624    506791  7686   7392   
minimap2-ensemble        3207988  4418   3132    503588  4543   10595  
'''

text = '''
strelka2-SNP             3206949  4038   4171    3202865  3342   8255   
vardict-SNP              3195487  27840  15633   3194962  26981  16158  
gatk-haplotype-SNP       3209073  8345   2047    3207496  7050   3624   
ensemble-SNP             3209551  4699   1569    3207988  4418   3132   
strelka2-indel        504322  30400  9861    503762  28823  10421  
vardict-indel         473815  16530  40368   475068  16338  39115  
gatk-haplotype-indel  506816  7771   7367    506791  7686   7392   
ensemble-indel        503714  4515   10469   503588  4543   10595  
'''

d = OrderedDict()

def _f_measure(b, prec, recall):
    assert (b**2 * prec + recall) != 0, (b, prec, recall)
    return (1 + b**2) * prec * recall / (b**2 * prec + recall)

for l in text.split('\n'):
    if l.strip():
        tokens = l.split()
        sn, stp, sfp, sfn, itp, ifp, ifn = l.split()

        snp = dict()
        ind = dict()
        
        snp['TP'] = int(stp)
        ind['TP'] = int(itp)
        snp['FN'] = int(sfn)
        ind['FN'] = int(ifn)
        snp['FP'] = int(sfp)
        ind['FP'] = int(ifp)
        
        snp['P']      = snp['TP'] + snp['FP']
        snp['Prec']   = float(snp['TP']) / snp['P'] if snp['P'] else 0
        snp['T']      = snp['TP'] + snp['FN']
        snp['Recall'] = float(snp['TP']) / snp['T'] if snp['T'] else 0

        ind['P']      = ind['TP'] + ind['FP']
        ind['Prec']   = float(ind['TP']) / ind['P'] if ind['P'] else 0
        ind['T']      = ind['TP'] + ind['FN']
        ind['Recall'] = float(ind['TP']) / ind['T'] if ind['T'] else 0

        snp['F1'] = _f_measure(1, snp['Prec'], snp['Recall'])
        snp['F2'] = _f_measure(2, snp['Prec'], snp['Recall'])
        ind['F1'] = _f_measure(1, ind['Prec'], ind['Recall'])
        ind['F2'] = _f_measure(2, ind['Prec'], ind['Recall'])

        d[sn] = (snp, ind)

import json
print(json.dumps(d, indent=4))



