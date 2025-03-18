import pandas as pd
geneinfo = pd.read_table('LD_snps_geneinfo.tsv')
LD = pd.read_csv('sig_mdtr.csv')
df = geneinfo.merge(LD, on = 'snp', how = 'left')
df.to_csv('sig_LD_geneinfo.tsv', sep = '\t', index = False)
