# start
# simulate PE
import glob
import os
from Bio import SeqIO
import statistics
import argparse
from scipy.stats import poisson
from statistics import mean
import pandas as pd
import numpy as np
import random
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
# set up path
snp_folder = args.i + '/removerec/'
summary_folder = args.i + '/summary/'
assembly_folder = args.i + '/co-assembly/'
allgenemut_file = glob.glob('%s/*.norecom.gooddepth.txt'%(snp_folder))
clonal_file = ('%s/../mergevcf/clonal_genelength_new.txt'%(snp_folder))
genus_num_cutoff = 3
simulation_round = 1000
pvalue_max = 0.05
max_depth_fold = 2  #fold change depth/gene copy for PE genes, max 2 copy of genes

try:
    os.mkdir(summary_folder)
except IOError:
    pass

################################################### Function ########################################################
def load_genemut(genemut_file,lineage_SNP):
    lineage = os.path.basename(genemut_file).split('.norecom')[0]
    species = lineage.split('_')[0]
    species_lineage.setdefault(species, set())
    species_lineage[species].add(lineage)
    lineage_SNP.setdefault(species, [0, dict()]) # No. SNPs, gene dict
    for lines in open(genemut_file,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\t')
            if lines_set[6] != 'None':
                # only SNPs on genes
                genename = lines_set[6]
                if genename in Cluster:
                    cluster = Cluster[genename]
                    # No. SNPs per gene cluster
                    lineage_SNP[species][0] += 1
                    lineage_SNP[species][1].setdefault(cluster, [0,set(),0,0])  # No. SNPs + genome set with SNPs, No. N SNPs, No. truncation SNPs
                    lineage_SNP[species][1][cluster][0] += 1
                    lineage_SNP[species][1][cluster][1].add(lineage)
                    if lines_set[8] in ['N', 'NN']:
                        # N SNPs
                        lineage_SNP[species][1][cluster][-1] += 1
                        if '*' in lines_set[9]:
                            # truncation
                            lineage_SNP[species][1][cluster][-2] += 1
    return lineage_SNP

def load_clonal(clonal_file):
    clonal = dict()
    for lines in open(clonal_file,'r'):
        if not lines.startswith('cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            species = lines_set[0].split('_')[0]
            ORFlength = int(lines_set[2])
            clonal.setdefault(species,[])
            clonal[species].append(ORFlength)
    return clonal

def pvalue_mutgene(SNP, ORFlength,gene_set):
    SNP_genome_set_all = set()
    pvalueset = set()
    mut_rate = float(SNP)/ORFlength
    num_lineage_in_species = len(species_lineage[species])
    allgenepvalue = ['gene\tpvalue\n']
    for cluster in Cluster_length:
        if cluster in gene_set:
            this_gene_length = np.median(Cluster_length[cluster])
            if cluster not in Cluster_multiplecopy_number:
                # cluster has max 1 gene in each lineage
                SNP_gene,lineage_num,truncation_gene,N_SNP_gene = gene_set.get(cluster,[0,set(),[1], 0,0]) # how many mutations
                num_strains_with_mut = (len(Cluster_count[cluster][species]) + 1)/num_lineage_in_species # num unique genes of the same cluster in a species (plus the lineage that was treated as blastn database)/ number of lineages
                # considering the prevalence of this gene among lineages
                pvalue = 1 - poisson.cdf(SNP_gene - 1,mut_rate*this_gene_length*num_strains_with_mut)# greater than and equal to
                allsum_details.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.5f\t%s\t%s\n'%(species,cluster,pvalue,SNP_gene,N_SNP_gene,truncation_gene,this_gene_length,num_strains_with_mut,mut_rate*1000,len(lineage_num)))
                allgenepvalue.append('%s\t%s\n' % (cluster, pvalue))
                if SNP_gene > 1 and len(lineage_num) > 1 and pvalue <= pvalue_max:
                    # potential PE
                    pvalueset.add(pvalue)
                    SNP_genome_set_all.add(len(lineage_num))
        else:
            allgenepvalue.append('%s\t%s\n' % (cluster, 1))
    pvalueset = list(pvalueset)
    pvalueset.sort(reverse=True)
    SNP_genome_set_all = list(SNP_genome_set_all)
    SNP_genome_set_all.sort()
    if len(pvalueset) > 0:
        foutput = open('%s/allgenes.poisson.pvalue.across.%s.details.txt' % (summary_folder, species), 'w')
        foutput.write(''.join(allgenepvalue))
        foutput.close()
    return [pvalueset, SNP_genome_set_all]

def mutation_sim(SNP, ORFlength,pvalueset,SNP_genome_set_all):
    cluster_num = []
    cluster_length_set = []
    mut_rate = float(SNP) / ORFlength
    num_lineage_in_species = len(species_lineage[species])
    allgenepvalue = ['simulation\tpvalue\n']
    print(len(Cluster_species[species]))
    for cluster in Cluster_species[species]:
        cluster_num.append(cluster)
        cluster_length_set.append(np.median(Cluster_length[cluster]))
    for i in range(0, simulation_round):
        cluster_mut = random.choices(cluster_num, weights=cluster_length_set, k=SNP)
        allclusters_mut = set(cluster_mut)
        allsim_cluster = dict()
        for cluster in cluster_num:
            if cluster in allclusters_mut and cluster not in Cluster_multiplecopy_number:
                # cluster has max 1 gene in each lineage
                SNP_cluster = cluster_mut.count(cluster)
                this_cluster_length = np.median(Cluster_length[cluster])
                num_strains_with_mut = (len(Cluster_count[cluster][species]) + 1) / num_lineage_in_species
                # considering the prevalence of this gene among lineages
                pvalue = 1 - poisson.cdf(SNP_cluster - 1, mut_rate * this_cluster_length * num_strains_with_mut) # greater than and equal to
                allgenepvalue.append('%s\t%s\n' % (i, pvalue))
                if SNP_cluster > 1:
                    # FP PE
                    for SNP_cluster_sub in range(2, SNP_cluster):
                        allsim_cluster.setdefault(SNP_cluster_sub, [])
                        allsim_cluster[SNP_cluster_sub].append(pvalue)
            else:
                allgenepvalue.append('%s\t%s\n' % (i, 1))
        for SNP_gene_sub in SNP_genome_set_all:
            for pvalue in pvalueset:
                num_genes_pass_pvalue = len([x for x in allsim_cluster.get(SNP_gene_sub,[]) if x <= pvalue])
                allsum.append('%s\t%s\t%s\t%s\t%s\n'%(species,SNP_gene_sub,pvalue,i,num_genes_pass_pvalue))
    if len(pvalueset) > 0:
        foutput = open('%s/allgenes.poisson.pvalue.across.simulation.%s.details.txt' % (summary_folder, species), 'w')
        foutput.write(''.join(allgenepvalue))
        foutput.close()
    return allsum

def load_gene_copy():
    allgenecopy = pd.DataFrame()
    allgenecopyfiles = glob.glob('%s/*.copynum.txt'%(assembly_folder))
    for genecopyfile in allgenecopyfiles:
        genecopy = pd.read_csv(genecopyfile, sep='\t')
        allgenecopy = allgenecopy.append(genecopy)
    # only considering single copy genes
    allgenecopy = allgenecopy[allgenecopy['copy_number']<max_depth_fold]
    return allgenecopy

def load_blastn(allblastnfiles,allgenecopy):
    Cluster = dict() # gene1 -> gene2
    Cluster_length = dict() # gene 2 -> length
    Cluster_count = dict() # gene 2 -> No. lineages with gene 2
    Cluster_species = dict() # species -> all gene 2
    Cluster_multiplecopy_number = [] # gene 2 with > 1 homologous genes in any lineage
    allgenename = list(allgenecopy['gene_name'])
    for blastnfile in allblastnfiles:
        species = os.path.basename(blastnfile).split('_')[0]
        lineage_name = os.path.basename(blastnfile).split('_%s'%(species))[0]
        for lines in open(blastnfile,'r'):
            lines_set = lines.split('\t')
            Gene1, Gene2, Identity, Length = lines_set[0:4]
            if Gene2 in allgenename and Gene1 in allgenename:
                # gene 2 single copy and gene 1 single copy
                Cluster_count.setdefault(Gene2, dict())
                Cluster_count[Gene2].setdefault(species, set())
                if lineage_name in Cluster_count[Gene2][species]:
                    # > 1 genes in lineage_name that belongs to gene 2 cluster
                    # should be removed
                    Cluster_multiplecopy_number.append(Gene2)
                else:
                    Length = int(Length)
                    Cluster.setdefault(Gene1, Gene2)
                    Cluster.setdefault(Gene2, Gene2)
                    Cluster_length.setdefault(Gene2, [Length])
                    Cluster_length[Gene2].append(Length)
                    Cluster_count[Gene2][species].add(lineage_name)
                    Cluster_species.setdefault(species, set())
                    Cluster_species[species].add(Gene2)
    return [Cluster, Cluster_length, Cluster_count, Cluster_species,Cluster_multiplecopy_number]

def compute_allcluster_species_length(species,Cluster_species,Cluster_multiplecopy_number,Cluster_length):
    ORFlength = 0
    for cluster in Cluster_species[species]:
        if cluster not in Cluster_multiplecopy_number:
            ORFlength += np.median(Cluster_length[cluster])
    return ORFlength

################################################### Main ########################################################
# load gene copy number
allgenecopy = load_gene_copy()

# load SNPs and genes with mutations
lineage_SNP_file = dict()
species_lineage = dict()
for genemut_file in allgenemut_file:
    lineage = os.path.basename(genemut_file).split('.norecom')[0]
    species = lineage.split('_')[0]
    lineage_SNP_file.setdefault(species,[])
    lineage_SNP_file[species].append(genemut_file)

# load non ORF length
clonal = load_clonal(clonal_file)
# simulation
allsum_details = []
allsum_details.append('species\tgene_name\tpvalue\tSNP_gene\tN_gene\ttruncation_gene\tgene_length\tprevalence\tmut_rate_1kbp\tlineage_set\n')
total_gene_length = 0
for species in lineage_SNP_file:
    genemut_filespecies = lineage_SNP_file[species]
    if len(genemut_filespecies) > 1:
        # load blastn results
        print('load blastn files for %s'%(species))
        allblastnfiles = glob.glob('%s/%s*clustercluster*clustercluster*.txt' % (assembly_folder, species))
        Cluster, Cluster_length, Cluster_count, Cluster_species, Cluster_multiplecopy_number = load_blastn(
            allblastnfiles,
            allgenecopy)
        # at least 2 lineages
        lineage_SNP = dict()
        for genemut_file in genemut_filespecies:
            print('processing %s' % (genemut_file))
            lineage_SNP = load_genemut(genemut_file, lineage_SNP)
        SNP, gene_set = lineage_SNP[species]
        print(SNP,species)
        if SNP > 2: #SNPs on gene > 1
                ORFlength = statistics.mean(clonal[species])
                print('original ORFlength %s for species %s'%(species,ORFlength))
                ORFlength = compute_allcluster_species_length(species, Cluster_species, Cluster_multiplecopy_number, Cluster_length)
                print('corrected ORFlength %s for species %s' % (species, ORFlength))
                pvalueset,SNP_genome_set_all = pvalue_mutgene(SNP, ORFlength, gene_set)
                print('finish compute poisson pvalue for species %s with ORFlength %s and SNPs %s PE pvalue set %s SNP genome set %s' % (species,ORFlength,SNP,pvalueset,SNP_genome_set_all))
                allsum = []
                allsum.append('species\tSNP_cutoff\tpvalue_cutoff\tsimulation_round\tFP\n')
                if len(pvalueset) > 0:
                # simulation
                    mutation_sim(SNP, ORFlength, pvalueset,SNP_genome_set_all)
                print('finish simulation %s' % (species))
                foutput = open('%s/allgenes.poisson.pvalue.%s.across.simulation.txt' % (summary_folder,species), 'w')
                foutput.write(''.join(allsum))
                foutput.close()

foutput = open('%s/allgenes.poisson.pvalue.across.txt'%(summary_folder), 'w')
foutput.write(''.join(allsum_details))
foutput.close()

################################################### END ########################################################
