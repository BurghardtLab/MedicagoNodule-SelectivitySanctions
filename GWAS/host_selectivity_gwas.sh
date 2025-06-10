#!/bin/bash 
#SBATCH --nodes=1 
#SBATCH --ntasks=20 
#SBATCH --mem=64GB 
#SBATCH --time=8:00:00 
#SBATCH --partition=open
 
salloc -N 1 -n 20 --mem=64GB -t 8:00:00 # use for interactive session

conda activate gwas_vcf # requires bcftools, vcftools, plink, plink2, gemma, bedtools, bedops, python, R

GENOTYPES="/storage/home/jus1990/burghardt/Jeremy/gwas2/requirements/medtr.A17.gnm5.div.Epstein_Burghardt_2022.Mt5_qual30_primitives_imputed.2021-09-29.vcf.gz" #SNPs file obatined from "wget https://data.legumeinfo.org/Medicago/truncatula/diversity/A17.gnm5.div.Epstein_Burghardt_2022/medtr.A17.gnm5.div.Epstein_Burghardt_2022.Mt5_qual30_primitives.2021-09-29.vcf.gz"
STRAINS="/storage/home/jus1990/burghardt/Jeremy/gwas2/requirements/strain_list.txt" #list of strains in phenotypes (i.e. column 1) # removed HM263 and HM272, not in Genotype file
PHENOTYPES="/storage/home/jus1990/burghardt/Jeremy/gwas2/requirements/MedicagoGWAS_traits_means.tsv" # TSV where first column = genotype list, phenotypes begin at column 2, no header
COVARIATES="/storage/home/jus1990/burghardt/Jeremy/gwas2/requirements/covariates.tsv" # TSV where first column = all 1's (intercept), covariates begin at column 2, no header
GENEFILE="/storage/home/jus1990/burghardt/Jeremy/gwas2/requirements/merged_genomic_CDS_INRA.gff" # https://medicago.toulouse.inra.fr/MtrunA17r5.0-ANR/downloads/1.9/20220708_MtrunA17r5.0-ANR-EGN-r1.9_Annotation.zip

bcftools view --samples-file "$STRAINS" "$GENOTYPES" --force-samples -o filtered_imputed.vcf # subset genotypes by samples, force samples if strain list has more than available in genotype file

bcftools view --max-alleles 2  -q 0.05:minor filtered_imputed.vcf -o filtered_biallelic_MAF005.vcf # Keep MAF 0.05 and biallelic positions only

sed -i 's/medtr.gnm5.A17.MtrunA17//g' filtered_biallelic_MAF005.vcf  #remove weird parts of chromosome names which creates errors in bed conversion downstream

#### index for downstream filtering ####
rm -f filtered_biallelic_MAF005.vcf.gz
bgzip filtered_biallelic_MAF005.vcf && tabix -p vcf filtered_biallelic_MAF005.vcf.gz

bcftools view --regions Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8 -o subset_chr1_8 filtered_biallelic_MAF005.vcf.gz # remove extra contigs from filtered vcf file.

# QC
rm -f QC_results
echo "No. of SNPs in filtered_imputed.vcf" > QC_results
bcftools query -f '%POS\n' filtered_imputed.vcf | wc -l >> QC_results # 40,907,386
echo "No. of SNPs in filtered_biallelic_MAF005.vcf.gz" >> QC_results
bcftools query -f '%POS\n' filtered_biallelic_MAF005.vcf.gz | wc -l >> QC_results # 8,514,373
echo "No. of SNPs in subset_chr1_8" >> QC_results
bcftools query -f '%POS\n' subset_chr1_8 | wc -l >> QC_results # 8,476,914
echo "No. of CHROMs in filtered_biallelic_MAF005.vcf.gz" >> QC_results
bcftools query -f '%CHROM\n' filtered_biallelic_MAF005.vcf.gz | uniq | wc -l >> QC_results # 40
echo "No. of CHROMs in subset_chr1_8" >> QC_results
bcftools query -f '%CHROM\n' subset_chr1_8 | uniq | wc -l >> QC_results # 8

# subset phenotype file
grep -f $STRAINS $PHENOTYPES > phenos_subset.fam  
# cut -f3 phenos_subset.fam > pred_prop_lobed.fam # unused
cut -f1,2,4 phenos_subset.fam > host_selectivity.fam # extract FID, IID, and phenotype from phenos_subset.fam file
# cut -f5 phenos_subset.fam > vessel_conductance.fam # unused

plink2 --vcf subset_chr1_8 -pheno host_selectivity.fam --double-id --no-parents --no-sex --no-pheno --make-bed --out "mdtr" # convert vcf to the Plink binary ped file format recognized by gemma. A direct vcf file will not work.

gemma -bfile mdtr -km 2 -gk 1 -o medtr_ksmat # Calculate k matrix. "-gk 1" for centered matrix, “-km 2” to accompany PLINK binary ped format

gemma -bfile mdtr -k ./output/medtr_ksmat.cXX.txt -lmm 4 -o medtr_ulmm_host_selectivity # gwas

# plot snps before annotation (see below) 

#### PLOTTING ###############################################################################################################################
# Use python to extract relevant data for plotting. R has a hard time read the large dataset. I trimmed the assoc.txt file to four columns required for plotting.

# import pandas as pd

# df = pd.read_table('output/medtr_ulmm_host_selectivity.assoc.txt') # change to relavent file if needed
# df = df[['rs', 'chr', 'ps', 'p_lrt']]
# df = df.rename(columns={"rs": "SNP", "chr" : "Chromosome", "ps" : "Position", "p_lrt" : "host_selectivity"}, errors="raise")
# df
# df.to_csv('output/medtr_ulmm_host_selectivity.csv', index = False) # change to relevant file if needed


# Use R to plot pvals
# library(CMplot)
# library(plyr)
# library(vroom)

# setwd('/storage/home/jus1990/burghardt/Jeremy/gwas2/imputed/output/') # change to relavent directory

# data <- vroom('medtr_ulmm_host_selectivity.csv') # change to relavent file
# head(data)

# SNPs <- c(data$SNP[data$host_selectivity<0.5e-7]) #set significance threshold (I aimed for ~10 snps)
# SNPs # 12 SNPs

# write(SNPs, file = "sig_snps.txt") ###### Used in annotation step ######

# CMplot(data, plot.type="m", multracks=FALSE, threshold=c(0.5e-7), threshold.lty=c(1), 
#        threshold.lwd=c(2), threshold.col=c("black","grey"), amplify=TRUE, 
#        signal.col= c("red"), signal.cex=1, file="jpg", file.name="final",
#        dpi=300, file.output=TRUE, verbose=TRUE, 
#        highlight=SNPs, highlight.text=SNPs, highlight.text.cex=1)

# CMplot(data,plot.type="q",box=FALSE,file="jpg",file.name="final",dpi=300,
#          conf.int=TRUE,conf.int.col=NULL,threshold.col="red",threshold.lty=2,
#          file.output=TRUE,verbose=TRUE,width=5,height=5)


#### ANNOTATION #############################################################################################################################

plink -bfile mdtr --r2 --ld-window 999999 --ld-window-kb 10000 --ld-window-r2 0.8 --ld-snp-list output/sig_snps.txt --out sig_mdtr # calculate snps withing LD (R2=0.8) of significant snps
sed -i 's/[ ][ ]*/ /g' sig_mdtr.ld # remove double tabs from LD output
cut -d ' ' -f 7 sig_mdtr.ld > LD_snps # extract snp list from column 7

gff2bed --max-mem 50G < $GENEFILE > genefile.bed # convert genefile to bed format. Tab delimited: chr, start, end, name, score, strand, attribute

vcftools --vcf subset_chr1_8 --snps LD_snps --recode --recode-INFO-all --out LD_snp_subset # Subset vcf to just 'LD snps' 
vcf2bed --max-mem 50G < LD_snp_subset.recode.vcf > LD_snp_subset.bed # convert LD subset vcf to bed format


range=10000
awk -v v=$range '{$2 = $2 - v; print}' LD_snp_subset.bed > LD_snp_subset_1.bed  # substract range from start position
awk -v v=$range '{$3 = $3 + v; print}' LD_snp_subset_1.bed > LD_snp_subset_2.bed # add range to end position
sed -i 's/ /\t/g' LD_snp_subset_2.bed

head LD_snp_subset.bed
head LD_snp_subset_2.bed

# Intersect bed files
bedtools intersect -a LD_snp_subset_2.bed -b genefile.bed -wb > LD_snps_geneinfo_1.bed # -wb Write the original entry in B for each overlap. Useful for knowing what A overlaps. Restricted by -f and -r.
bedtools intersect -a LD_snp_subset_2.bed -b genefile.bed -v > LD_snps_geneinfo_2.bed # -v 	Only report those entries in A that have no overlap in B. Restricted by -f and -r.
rm -f LD_snps_geneinfo.bed
cat LD_snps_geneinfo*.bed > LD_snps_geneinfo.bed # merged intersects

sed 's/;GO/,GO/g' LD_snps_geneinfo.bed > LD_snps_geneinfo.tsv # convert ill placed ';'
sed -i 's/;/\t/g' LD_snps_geneinfo.tsv # expand 'attirbute' column in bed file to multiple columns

wc -l output/sig_snps.txt # 12 Significant snps (0.5e-7)
wc -l LD_snps # 309 LD_snps
wc -l LD_snps_geneinfo.bed # 1432 LD_snps_geneinfo.bed

# Transfer outputs to local computer
# scp -r jus1990@submit.hpc.psu.edu://storage/home/jus1990/burghardt/Jeremy/gwas2/imputed/LD_snps_geneinfo.tsv  .
# scp -r jus1990@submit.hpc.psu.edu://storage/home/jus1990/burghardt/Jeremy/gwas2/imputed/sig_mdtr.ld .
# scp -r jus1990@submit.hpc.psu.edu://storage/home/jus1990/burghardt/Jeremy/gwas2/imputed/output/*jpg .
# scp -r jus1990@submit.hpc.psu.edu://storage/home/jus1990/burghardt/Jeremy/gwas2/imputed/output/medtr_ulmm_host_selectivity* .


#Alternate GENEFILE
# GENEFILE="/storage/home/jus1990/burghardt/Jeremy/gwas2/requirements/genomic_CDS.gff" # grep CDS from genomic.gff, https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_003473485.1/download?include_annotation_type=GENOME_GFF
