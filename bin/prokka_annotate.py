# start
import os,glob
from Bio import SeqIO
import argparse
from scipy.stats import poisson
from statistics import mean
import pandas as pd
import numpy as np
import random
from scipy.stats import binom

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="a folder to store all output",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round3/',
                      metavar='output/')

################################################## Definition ########################################################
args = parser.parse_args()
workingdir=os.path.abspath(os.path.dirname(__file__))

################################################### Function ########################################################
def prokka_summary(prokka_file):
    # remove useless information
    try:
        f1 = open(prokka_file + '.brief','r')
    except IOError:
        alloutput = []
        for lines in open(prokka_file, 'r'):
            if lines.startswith('##FASTA'):
                break
            if not lines.startswith('#'):
                alloutput.append(lines)
        f1 = open(prokka_file + '.brief','w')
        f1.write(''.join(alloutput))
        f1.close()

def load_prokka(prokka_file):
    # load all genes prokka annotation
    lineage = os.path.basename(prokka_file).split('.noHM.fna')[0]
    species = lineage.split('_')[0]
    prokka = pd.read_csv(prokka_file + '.brief', sep='\t', header=None)
    prokka['product'] = ''
    prokka['gene'] = ''
    for i in prokka.index:
        if ';product=' in str(prokka.loc[i, 8]):
            prokka.loc[i, 'product'] = prokka.loc[i, 8].split(';product=')[1]
        if 'gene=' in str(prokka.loc[i, 8]):
            prokka.loc[i, 'gene'] = prokka.loc[i, 8].split('gene=')[1].split(';')[0]
    prokka = prokka[~prokka['product'].isin(['hypothetical protein', ''])]
    prokka = prokka.loc[:, [0, 'product', 'gene']]
    prokka['lineage_gene'] = ['%s__%s'%(lineage,x) for x in prokka[0]]
    prokka['species_gene'] = ['%s__%s' % (species, x) for x in prokka[0]]
    return prokka

def load_mutation_genes(gene_within_file,gene_across_file,prokka_all):
    gene_within = pd.read_csv(gene_within_file, sep='\t')
    gene_within['lineage_gene'] = gene_within['lineage'] + '__' + gene_within['gene_name']
    gene_across = pd.read_csv(gene_across_file, sep='\t')
    gene_across['species_gene'] = gene_across['species'] + '__' + gene_across['gene_name']
    print(gene_within.head())
    print(prokka_all.head())
    gene_within = gene_within.merge(prokka_all, left_on='lineage_gene', right_on='lineage_gene', how='left')
    gene_within.to_csv(gene_within_file.replace('.txt','_prokka.txt'), sep='\t',
                              index=False)
    gene_across = gene_across.merge(prokka_all, left_on='species_gene', right_on='species_gene', how='left')
    gene_across.to_csv(gene_across_file.replace('.txt', '_prokka.txt'), sep='\t',
                       index=False)

################################################### Main ########################################################
gene_within_file = glob.glob('%s/summary/*within*cutoff.txt'%(args.i))[0]
gene_across_file = glob.glob('%s/summary/*across*cutoff.txt'%(args.i))[0]
# load all prokka
allprokkafiles = glob.glob('%s/co-assembly/*.gff'%(args.i))
prokka_all = pd.DataFrame()
for prokka_file in allprokkafiles:
    print('load prokka for %s'%(prokka_file))
    prokka_summary(prokka_file)
    prokka_all = prokka_all.append(load_prokka(prokka_file))

# add prokka annotation to genes with mutations
print('add prokka annotation to genes with mutations')
load_mutation_genes(gene_within_file,gene_across_file,prokka_all)
################################################### END ########################################################
