# MedicagoNodule-SelectivitySanctions

This repository contains the majority of the data and code presented within the Burghardt et al., 2025 manuscript.

R Folder (contains analysis of plant phenotypes and completed in R Studio)
- Data (input data for analyses)
    -epstein_et_al_2022 (data from HapMap Select & Resequence)
        - epstein_et_al_2022_S9.csv - Relative strain frequencies from HapMap Select & Resequence nodule pools
    -HapMap-SpanFran2-Spring2020
        - df_summary.csv - summarized image analysis data of nodule pools from HapMap Select & Resequence nodule pools
        - S&RGWAS_PlantMetaData&Nodules_WInter2020.csv - experimental metadata for HapMap Select & Resequence grow
    - freq_C68.txt - 
    - NodulePools_SizeData.txt - 
    - NoduleSize_HarvestData.txt - 
    - NoduleSizeClassInfo.txt - 
    - SingleNoduleColonies_11Apr18.txt - 
    - SingleStrain_phenotype_summary.tsv - 
    
- Figures (figures presented within manuscript)
    - Names of PDFs align with manuscript figure captions
- Source (R scripts containing analyses)
    - NodSize-Rcode.R
    - NoduleMorphology_HapMap.R - Analyses of HapMap nodule morphologies, Host Selectivity simulations, and relationships between nodule morphologies and Host Selectivity

GWAS Folder (contains script and output files of the presented Genome-Wide Association Study)
- scripts (shell and R scripts for executing GWAS and visualizing results)
    - host_selectivity_gwas.sh
    - manhattan_plot.R
    - merge_sig_LD_geneinfo.py
    - process_GENEFILE.sh
- output (GWAS results and plots)
    - host_selectivity_QQplot_final.jpg
    - host_selectivity_Rect_Manhtn_final.jpg
    - LD_snps_geneinfo.tsv
    - medtr_ulmm_host_selectivity.assoc.txt
    - medtr_ulmm_host_selectivity.csv
    - medtr_ulmm_host_selectivity.log.txt
    - sig_LD_geneinfo_pval_beta_compiled.xlsx
    - sig_LD_geneinfo_pval_beta.tsv
    - sig_LD_geneinfo.tsv
    - sig_mdtr.csv
    - sig_mdtr.ld


