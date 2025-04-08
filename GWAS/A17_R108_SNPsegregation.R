library(tidyverse)

df <- read_tsv("sig_LD_geneinfo_pval_beta.tsv")

sig.list <- c("CM010654.1-26952852-C", "CM010650.1-51530113-A", "CM010651.1-51084261-T", 
"CM010653.1-25128406-A", "CM010651.1-51080330-T", "CM010654.1-23572058-T",
"CM010653.1-32706149-T", "CM010650.1-45833295-C", "CM010649.1-16046163-T", 
"CM010654.1-23572805-T")

sig.list <- paste0("freebayes-var-", sig.list)

# Filter to SNPs of interest and HM101 (A17) and HM340 (R108)
df <- df %>% filter(snp %in% sig.list) %>% 
  select(snp, pvalue, beta, HM101) %>% unique()

df <- read.csv("HapMap_SNPs.csv")

df <- df %>% filter(X %in% sig.list) %>% 
  select(X, HM101, HM340)
