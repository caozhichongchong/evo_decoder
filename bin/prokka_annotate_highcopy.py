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
    #prokka['species_gene'] = ['%s__%s' % (species, x) for x in prokka[0]]
    return prokka

def load_mutation_genes(gene_copy_files,prokka_all):
    allgenes = pd.DataFrame()
    for files in gene_copy_files:
        lineage_name = os.path.basename(files).split('.noHM.fna')[0]
        allgenessub = pd.read_csv(files, sep='\t')
        allgenessub = allgenessub[allgenessub['copy_number']>=2]
        allgenessub['lineage'] = lineage_name
        allgenes = allgenes.append(allgenessub)
    allgenes['lineage_gene'] = allgenes['lineage'] + '__' + allgenes['gene_name']
    print(allgenes.head())
    print(prokka_all.head())
    allgenes = allgenes.merge(prokka_all, left_on='lineage_gene', right_on='lineage_gene', how='left')
    allgenes.to_csv('%s/summary/allgenes.highcopynum.prokka.txt'%(args.i), sep='\t',
                              index=False)

################################################### Main ########################################################
gene_copy_files = glob.glob('%s/co-assembly/BaFr*.copynum.txt'%(args.i))

# load all prokka
allprokkafiles = glob.glob('%s/co-assembly/BaFr*.gff'%(args.i))
prokka_all = pd.DataFrame()
for prokka_file in allprokkafiles:
    print('load prokka for %s'%(prokka_file))
    prokka_summary(prokka_file)
    prokka_all = prokka_all.append(load_prokka(prokka_file))

# add prokka annotation to genes with mutations
print('add prokka annotation to genes with mutations')
load_mutation_genes(gene_copy_files,prokka_all)
################################################### END ########################################################
