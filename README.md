### RhizobiaSanctioning-NoduleSize

Repository for Burghardt et al., 2025 data and associated analyses including 
plant phenotypes and summaries of GWAS outputs

# Source

NoduleMorphology_HapMap.R - Contains code for Medicago HapMap panel analysis with host selectivity simulations, analysis of nodule traits, and relationships between nodule traits and host selectivity producing Figure 2 and Table S1.

NodSize_A17_R108.R - Contains code for comparisons between A17 and R108 within Figures 3, 4, S3, S4, and S5 and associated analyses producing Tables S4, S5, and S6. 

# Data

NodulePools_SizeData.txt -  Overall characteristics of A17 and R108 nodule size pools (e.g. number of nodules and wet weight)

SingleNoduleColonies.txt - Colony forming units estimated from individual crushed A17 and R108 nodules of each size class. ~ 3 replicate dilution series were conducted.

NoduleSizeClassInfo.txt - Measurements A17 and R108 nodule length and number of branches from 10 randomly chosen nodules of each Sized pool

freq_C68.txt - Frequencies of 68 strains within nodules communnities of A17 and R108 small, medium, and large pools

SingleStrain_phenotype_summary.tsv -  Single-strain plant growth benefit data for 68 strains from a previously published experiment (Burghardt 2018)

epstein_et_al_2022 - Contains previously published data describing the frequencies of 88 strains within initial inocula and nodule pools of the HapMap panel  (HapMap Select & Resequence Experiment (Epstein et al., 2022))

HapMap-SpanFran2-Spring2020 - Contains nodule phenotypes recorded at time of harvest for the 2020 HapMap Select & Resequence experiment (S&RGWAS_PlantMetaData&Nodules_WInter2020.csv) as well as the summary dataset of ImageJ analysis output measures (df_summary.csv)

# Figures

Contains original figure pdf files from R Studio as well as final edited files in Adobe Illustator for manuscript figures

# GWAS

Includes summaries of host selectivity GWAS pipeline (host_selectivity_gwas.sh, medtr_ulmm_host_selectivity.log.txt), GWAS output, and output analysis (merge_sig_LD_geneinfo.py,manhattan_plot.R, A17_R108_SNPsegregation.R)

Significant SNP information is mainly within sig_LD_geneinfo_pval_beta.tsv

This folder does not contain larger GWAS output files such as

medtr_ulmm_host_selectivity.assoc.txt - SNPs associated with host selectivity

medtr_ulmm_host_selectivity.csv - ulmm output of host selectivity means associated with SNPs

