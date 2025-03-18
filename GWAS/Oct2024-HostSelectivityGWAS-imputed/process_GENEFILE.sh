## HOW to process merged_genomic_CDS_INRA.gtf
grep 'CDS' MtrunA17r5.0-ANR-EGN-r1.9.gtf > genomic_CDS_INRA.gff2
grep 'A17Chr' genomic_CDS_INRA.gff2 > genomic_CDS_INRA.gff
sed -i 's/,/|/g' genomic_CDS_INRA.gff # replace ,
sed -i 's/;/\t/g' genomic_CDS_INRA.gff # remove ;
sed -i 's/[ ][ ]*/ /g' genomic_CDS_INRA.gff # remove double spaces
sed -i 's/gene_biotype //g' genomic_CDS_INRA.gff # remove tags
sed -i 's/ gene_id //g' genomic_CDS_INRA.gff
sed -i 's/ gene_name //g' genomic_CDS_INRA.gff
sed -i 's/ locus_tag //g' genomic_CDS_INRA.gff
sed -i 's/ transcript_id /\t/g' genomic_CDS_INRA.gff
sed -i 's/ product //g' genomic_CDS_INRA.gff
sed -i 's/"//g' genomic_CDS_INRA.gff # remove quotes
echo -e "seqname\tsource\tfeature\tstart\tstop\tscore\tstrand\tframe\tgene_biotype\tgene_id\tgene_name\tlocus_tag\tproduct\ttranscript_id" | cat - genomic_CDS_INRA.gff > genomic_CDS_INRA.gff2 # add new header
awk 'BEGIN{FS=OFS="\t"} {NF=13} 1' genomic_CDS_INRA.gff2 > genomic_CDS_INRA.gff # fix tabs


#### Merge with summary tsv and create new attibute column
# # Python
# import pandas as pd 
# gff = pd.read_table('genomic_CDS_INRA.gff')
# summ = pd.read_table('MtrunA17r5.0-ANR-EGN-r1.9.summary.tsv')
# df = gff.merge(summ, on = 'locus_tag', how = 'left')
# df.to_csv('genomic_CDS_INRA.tsv', sep="\t") 
# df = pd.DataFrame(df)

# merged_column = pd.Series([';'.join(f'{col}={val}' for col, val in zip(df.columns[9:], row[9:])) for _, row in df.iterrows()]) # Merge columns with '=' as separator.
# df['attribute'] = merged_column

# df1 = df[['seqname', 'source', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'attribute']]
# df1['seqname'] = df['seqname'].str.replace('MtrunA17', '')
# df1.to_csv('merged_genomic_CDS_INRA.gff', sep='\t', index=False, header=False)
