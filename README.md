# evo_decoder
1) prep meta data\
genome  lineage fastq   fasta\
bq_BL_g0001     BL_clustercluster13.donor.bq    PATH_TO_1.fastq    PATH_TO_ASSEMBLY
2) run evo_decoder :)\
`python mapping_WGS.py -p PATH_GENOME_META_DATA -o OUTPUT_DIR -s SCRIPT_DIR`\
`python vcf_process.py -i OUTPUT_DIR -s SCRIPT_DIR`\
`sh SCRIPT_DIR/vcf_filter.sh`\
`python SNP_merge.py -i OUTPUT_DIR`\
`python cluster_length_new.py -i OUTPUT_DIR`\
NOTEBOOk PE_dnds.ipynb\
`python PE_gene_pvalue.py -i OUTPUT_DIR`\
`python PE_gene_pvalue_across.py -i OUTPUT_DIR`
3) extra prokka annotation\
upload *cutoff.txt to OUTPUT_DIR/summary/\
`python prokka_annotate.py`
 

