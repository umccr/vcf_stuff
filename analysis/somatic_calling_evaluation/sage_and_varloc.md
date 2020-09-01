```
BCBIO_TNAME=T_SRR7890902_20pc
TBAM=/g/data/gx8/projects/Saveliev_SEQCII/bcbio/samples/final/T_SRR7890902_20pc/T_SRR7890902_20pc-ready.bam
NBAM=/g/data/gx8/projects/Saveliev_SEQCII/bcbio/samples/final/N_SRR7890889/N_SRR7890889-ready.bam
TNAME=SEQCII_20pc
NNAME=SEQCII_N
```

# sage

```
cd /g/data/gx8/projects/Saveliev_SEQCII/sage

java -Xms4G -Xmx32G -cp sage-2.2-jar-with-dependencies.jar com.hartwig.hmftools.sage.SageApplication \
    -threads 16 \
    -reference $NNAME -reference_bam $NBAM \
    -tumor $TNAME -tumor_bam $TBAM \
    -assembly hg38 \
    -ref_genome /g/data/gx8/extras/umccrise/genomes/hg38/hg38.fa \
    -hotspots /g/data/gx8/extras/umccrise/genomes/hg38/hotspots/merged.vcf.gz \
    -panel_bed /g/data/gx8/extras/umccrise/genomes/hg38/hmf/coding_regions.bed \
    -high_confidence_bed /g/data/gx8/extras/umccrise/genomes/hg38/hmf/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed.gz \
    -out $TNAME.sage.vcf.gz

sage_vcf=$TNAME.sage.vcf.gz
annotated_vcf=$TNAME.sage.pon_ann.vcf.gz
pon_filtered_vcf=$TNAME.sage.pon_ann.pon_filt.vcf.gz
bcftools annotate -a SageGermlinePon.hg38.98x.vcf.gz -c PON_COUNT,PON_MAX $sage_vcf -O z -o $annotated_vcf
bcftools filter -e 'PON_COUNT!="." && INFO/TIER="HOTSPOT" && PON_MAX>=5 && PON_COUNT >= 2' -s PON -m+ $annotated_vcf -O u | \
	bcftools filter -e 'PON_COUNT!="." && INFO/TIER="PANEL" && PON_MAX>=5 && PON_COUNT >= 2' -s PON -m+ -O u | \
	bcftools filter -e 'PON_COUNT!="." && INFO/TIER!="HOTSPOT" && INFO/TIER!="PANEL" && PON_COUNT >= 5' -s PON -m+ -O z -o $pon_filtered_vcf
bcftools view -f.,PASS $TNAME.sage.vcf.gz -Oz -o $TNAME.sage.PASS.vcf.gz
bcftools view -f.,PASS $TNAME.sage.pon_ann.pon_filt.vcf.gz -Oz -o $TNAME.sage.pon_ann.pon_filt.PASS.vcf.gz

pcgr_prep $TNAME.sage.pon_ann.pon_filt.PASS.vcf.gz -tn $TNAME -nn $NNAME -o $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.vcf.gz
pcgr_prep $TNAME.sage.PASS.vcf.gz -tn $TNAME -nn $NNAME -o $TNAME.sage.PASS.pcgr_prep.vcf.gz
bcftools view -T /g/data3/gx8/extras/hg38_noalt_noBlacklist.bed $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.vcf.gz -Oz -o $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.vcf.gz
anno_somatic_vcf -tn $TNAME -nn $NNAME $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.vcf.gz -o $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.vcf.gz -g hg38 
filter_somatic_vcf -tn $TNAME -nn $NNAME $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.vcf.gz -o $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.filt.vcf.gz -g hg38
bcftools view -s $TNAME $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.filt.vcf.gz -Oz -o $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.filt.1sample.vcf.gz
bcftools view -f.,PASS $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.filt.1sample.vcf.gz -Oz -o $TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.filt.1sample.PASS.vcf.gz

```

# dragen and bcbio - SAGE PoN

```
# ensemble
cd /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio
sage_vcf=${TNAME}-ensemble-annotated.vcf.gz
annotated_vcf=${TNAME}.sage_pon.vcf.gz
pon_filtered_vcf=${TNAME}.sage_pon.filt.vcf.gz
pon_hard_filtered_vcf=${TNAME}.sage_pon.filt.pass.vcf.gz
bcftools annotate -a SageGermlinePon.hg38.98x.vcf.gz -c PON_COUNT,PON_MAX $sage_vcf -O z -o $annotated_vcf
bcftools filter $annotated_vcf -e 'PON_COUNT!="." && PON_MAX>=5 && PON_COUNT >= 5' -s PON -m+ -O z -o $pon_filtered_vcf
bcftools view -f.,PASS $pon_filtered_vcf -Oz -o $pon_hard_filtered_vcf

# strelka
cd /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio
bcftools view -f.,PASS ../../bcbio/samples/final/${BCBIO_TNAME}/${BCBIO_TNAME}-strelka2.vcf.gz -Oz -o ${TNAME}-strelka.vcf.gz
tabix ${TNAME}-strelka.vcf.gz
sage_vcf=${TNAME}-strelka.vcf.gz
annotated_vcf=${TNAME}-strelka.sage_pon.vcf.gz
pon_filtered_vcf=${TNAME}-strelka.sage_pon.filt.vcf.gz
pon_hard_filtered_vcf=${TNAME}-strelka.sage_pon.filt.pass.vcf.gz
bcftools annotate -a SageGermlinePon.hg38.98x.vcf.gz -c PON_COUNT,PON_MAX $sage_vcf -O z -o $annotated_vcf
bcftools filter $annotated_vcf -e 'PON_COUNT!="." && PON_MAX>=5 && PON_COUNT >= 5' -s PON -m+ -O z -o $pon_filtered_vcf
bcftools view -f.,PASS $pon_filtered_vcf -Oz -o $pon_hard_filtered_vcf

# dragen
cd /g/data/gx8/projects/Saveliev_SEQCII/validation/dragen
sage_vcf=${TNAME}.PREP.NOALT.vcf.gz
annotated_vcf=${TNAME}.sage_pon.vcf.gz
pon_filtered_vcf=${TNAME}.sage_pon.filt.vcf.gz
pon_hard_filtered_vcf=${TNAME}.sage_pon.filt.pass.vcf.gz
tabix $sage_vcf
bcftools annotate -a SageGermlinePon.hg38.98x.vcf.gz -c PON_COUNT,PON_MAX $sage_vcf -O z -o $annotated_vcf
bcftools filter $annotated_vcf -e 'PON_COUNT!="." && PON_MAX>=5 && PON_COUNT >= 5' -s PON -m+ -O z -o $pon_filtered_vcf
bcftools view -f.,PASS $pon_filtered_vcf -Oz -o $pon_hard_filtered_vcf
```



# varlociraptor

```
cd /g/data/gx8/projects/Saveliev_SEQCII/varlociraptor

varlociraptor preprocess variants /g/data/gx8/extras/umccrise/genomes/hg38/hg38.fa \
	--bam $TBAM \
	--candidates ../bcbio/samples/final/2019-11-25_samples/SEQCII_50pc-ensemble-annotated.vcf.gz \
	--output $TNAME_varloc.bcf

varlociraptor preprocess variants /g/data/gx8/extras/umccrise/genomes/hg38/hg38.fa \
	--bam $NBAM \
	--candidates ../bcbio/samples/final/2019-11-25_samples/SEQCII_50pc-ensemble-annotated.vcf.gz \
	--output $NNAME_varloc.bcf

varlociraptor call variants tumor-normal --purity 0.5 \
	--tumor $TNAME_varloc.bcf \
	--normal $NNAME_varloc.bcf \
	> $TNAME_varloc.calls.bcf

cat ${TNAME}_varloc.calls.bcf | varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds positive \
	| bcftools view -Oz -o ${TNAME}_varloc.calls.positive.vcf.gz
cat ${TNAME}_varloc.calls.bcf | varlociraptor filter-calls posterior-odds --events SOMATIC_TUMOR --odds strong \
	| bcftools view -Oz -o ${TNAME}_varloc.calls.strong.vcf.gz

varlociraptor filter-calls control-fdr ${TNAME}_varloc.calls.bcf --events SOMATIC_TUMOR --fdr 0.05 --var SNV \
	| bcftools view -Oz -o ${TNAME}_varloc.calls.fdr005.vcf.gz
varlociraptor filter-calls control-fdr ${TNAME}_varloc.calls.bcf --events SOMATIC_TUMOR --fdr 0.1 --var SNV \
	| bcftools view -Oz -o ${TNAME}_varloc.calls.fdr01.vcf.gz

```

# Evaluation

```
mkdir /g/data/gx8/projects/Saveliev_SEQCII/validation/eval_sage_varloc_50pc
cd /g/data/gx8/projects/Saveliev_SEQCII/validation/eval_sage_varloc_50pc
TNAME=SEQCII_50pc

# 1. Bcbio                          
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio/$TNAME-ensemble-annotated.vcf.gz bcbio.vcf.gz
# 3. Bcbio, umccrised              
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio/umccrised/$TNAME.PASS.vcf.gz bcbio.umccrised.vcf.gz
# 2. Dragen                         
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/dragen/$TNAME.PREP.NOALT.vcf.gz dragen.vcf.gz
# 4. Dragen, umccrised              
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/dragen/$TNAME.PREP.NOALT.ANNO.FILT.PASS.vcf.gz dragen.umccrised.vcf.gz
# 5. SAGE
ln -s /g/data/gx8/projects/Saveliev_SEQCII/sage/$TNAME.sage.PASS.vcf.gz sage.vcf.gz
# 6. SAGE + PON
ln -s /g/data/gx8/projects/Saveliev_SEQCII/sage/$TNAME.sage.pon_ann.pon_filt.PASS.vcf.gz sage.PoN.vcf.gz
# 7. SAGE, umccrised
# 8. SAGE + PON, umccrised
ln -s /g/data/gx8/projects/Saveliev_SEQCII/sage/$TNAME.sage.pon_ann.pon_filt.PASS.pcgr_prep.noalt.anno.filt.1sample.PASS.vcf.gz sage.PoN.umccrised.vcf.gz
# 9. Bcbio, varlociraptor positive 
ln -s /g/data/gx8/projects/Saveliev_SEQCII/varlociraptor/${TNAME}_varloc.calls.positive.vcf.gz bcbio.varlociraptor.positive.vcf.gz
# 10. Bcbio, varlociraptor strong   
ln -s /g/data/gx8/projects/Saveliev_SEQCII/varlociraptor/${TNAME}_varloc.calls.strong.vcf.gz bcbio.varlociraptor.strong.vcf.gz
# 11. bcbio. varloctiraptor FDR 0.05
ln -s /g/data/gx8/projects/Saveliev_SEQCII/varlociraptor/${TNAME}_varloc.calls.fdr005.vcf.gz bcbio.varlociraptor.fdr005.vcf.gz
# 12. bcbio. varloctiraptor FDR 0.1
ln -s /g/data/gx8/projects/Saveliev_SEQCII/varlociraptor/${TNAME}_varloc.calls.fdr01.vcf.gz  bcbio.varlociraptor.fdr01.vcf.gz

TRUTH=/g/data/gx8/projects/Saveliev_SEQCII/validation/truth/truth_merged.ANNO.HighConf.MAKE_PASS.vcf.gz
TRUTH_AF10=/g/data/gx8/projects/Saveliev_SEQCII/validation/truth/truth_merged.ANNO.HighConf.MAKE_PASS.AF10.vcf.gz
eval_vcf $TRUTH_AF10 bcbio*.vcf.gz dragen* sage* -o eval_AF10

             SNP                                  
                     TP     FP    FN Recall   Prec     F2
           bcbio  29294  86231   391 98.68% 25.36% 62.52%
 bcbio.umccrised  26507   8215  3178 89.29% 76.34% 86.36%
          dragen  29114  74549   571 98.08% 28.09% 65.45%
            sage  28286  80031  1399 95.29% 26.11% 62.29%
        sage.PoN  28012  11456  1673 94.36% 70.97% 88.53%

           INDEL
                     TP     FP    FN Recall   Prec     F2
           bcbio    837  13162    12 98.59%  5.98% 24.06%
 bcbio.umccrised    732   1275   117 86.22% 36.47% 67.74%
          dragen    841  17350     8 99.06%  4.62% 19.48%
            sage    794  13862    55 93.52%  5.42% 21.99%
        sage.PoN    763   1897    86 89.87% 28.68% 63.00%
                                    
```

```
                             SNP TP     FP    FN Recall   Prec     F2
                       bcbio  29294  86231   391 98.68% 25.36% 62.52%
  bcbio.varlociraptor.fdr005  29294  86224   391 98.68% 25.36% 62.53%
   bcbio.varlociraptor.fdr01  29294  86224   391 98.68% 25.36% 62.53%
bcbio.varlociraptor.positive  29179  84627   506 98.30% 25.64% 62.74%
  bcbio.varlociraptor.strong  29238  85323   447 98.49% 25.52% 62.66%
```

```
mkdir /g/data/gx8/projects/Saveliev_SEQCII/validation/eval_sage_varloc_50pc
cd /g/data/gx8/projects/Saveliev_SEQCII/validation/eval_sage_varloc_50pc

SEQCII, 50pc                                         
SNP                 TP     FP    FN Recall   Prec     F2
          bcbio  29294  86231   391 98.68% 25.36% 62.52%
        strelka  29192  86126   493 98.34% 25.31% 62.36%
         dragen  29114  74549   571 98.08% 28.09% 65.45%
           sage  28286  80031  1399 95.29% 26.11% 62.29%
  bcbio.sagePoN  29049  16478   636 97.86% 63.81% 88.42%
strelka.sagePoN  28952  14694   733 97.53% 66.33% 89.15%
 dragen.sagePoN  28875  14238   810 97.27% 66.98% 89.20%
   sage.sagePoN  28012  11456  1673 94.36% 70.97% 88.53%

INDEL               TP     FP    FN Recall   Prec     F2
          bcbio    837  13162    12 98.59%  5.98% 24.06%
        strelka    824  15413    25 97.06%  5.07% 20.99%
         dragen    841  17350     8 99.06%  4.62% 19.48%
           sage    794  13862    55 93.52%  5.42% 21.99%
  bcbio.sagePoN    811   2196    38 95.52% 26.97% 63.33%
strelka.sagePoN    796   2169    53 93.76% 26.85% 62.57%
 dragen.sagePoN    810   3526    39 95.41% 18.68% 52.38%
   sage.sagePoN    763   1897    86 89.87% 28.68% 63.00%
```

```
mkdir /g/data/gx8/projects/Saveliev_SEQCII/validation/eval_sage_varloc_20pc
cd /g/data/gx8/projects/Saveliev_SEQCII/validation/eval_sage_varloc_20pc

ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio/${TNAME}.sage_pon.filt.pass.vcf.gz          bcbio.sagePoN.vcf.gz            
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio/${TNAME}-ensemble-annotated.vcf.gz          bcbio.vcf.gz            
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/dragen/${TNAME}.sage_pon.filt.pass.vcf.gz         dragen.sagePoN.vcf.gz             
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/dragen/${TNAME}.PREP.NOALT.vcf.gz                 dragen.vcf.gz     
ln -s /g/data/gx8/projects/Saveliev_SEQCII/sage/${TNAME}.sage.pon_ann.pon_filt.PASS.pcgr_prep.vcf.gz    sage.PoN.vcf.gz                  
ln -s /g/data/gx8/projects/Saveliev_SEQCII/sage/${TNAME}.sage.PASS.pcgr_prep.vcf.gz                     sage.vcf.gz 
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio/${TNAME}-strelka.sage_pon.filt.pass.vcf.gz  strelka.sagePoN.vcf.gz                    
ln -s /g/data/gx8/projects/Saveliev_SEQCII/validation/bcbio/${TNAME}-strelka.vcf.gz                     strelka.vcf.gz 

eval_vcf $TRUTH strelka.vcf.gz strelka.sagePoN.vcf.gz bcbio.vcf.gz bcbio.sagePoN.vcf.gz dragen.vcf.gz dragen.sagePoN.vcf.gz sage.vcf.gz sage.PoN.vcf.gz -o eval_sage_pon

SEQCII, 20pc    
SNP                    TP     FP   FN Recall   Prec     F2
            bcbio   26522  78147 3163 89.34% 25.34% 59.36%
          strelka   25262  76189 4423 85.10% 24.90% 57.36%
           dragen   25361  66613 4324 85.43% 27.57% 60.18%
             sage   20133  62716 9552 67.82% 24.30% 49.94%
    bcbio.sagePoN   26282  10894 3403 88.54% 70.70% 84.28%
  strelka.sagePoN   25032   9195 4653 84.33% 73.14% 81.82%
   dragen.sagePoN   25135   9069 4550 84.67% 73.49% 82.17%
     sage.sagePoN   19894   6662 9791 67.02% 74.91% 68.46%
   
INDEL                  TP     FP   FN Recall   Prec     F2
            bcbio     733   9893  116 86.34%  6.90% 26.14%
          strelka     642  10898  207 75.62%  5.56% 21.49%
           dragen     725  12644  124 85.39%  5.42% 21.62%
             sage     508   7148  341 59.84%  6.64% 22.98%
    bcbio.sagePoN     715   1581  134 84.22% 31.14% 62.81%
  strelka.sagePoN     624   1388  225 73.50% 31.01% 57.69%
   dragen.sagePoN     701   2352  148 82.57% 22.96% 54.35%
     sage.sagePoN     491    948  358 57.83% 34.12% 50.78%
```

















