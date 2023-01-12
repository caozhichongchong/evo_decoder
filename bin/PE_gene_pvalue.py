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
# set up path
snp_folder = args.i + '/removerec/'
MW_folder = args.i + '/vcf/MV_cov/'
assembly_folder = args.i + '/co-assembly/'
summary_folder = args.i + '/summary/'
allgenemut_file = glob.glob('%s/*.norecom.gooddepth.txt'%(snp_folder))
clonal_file = ('%s/../mergevcf/clonal_genelength_new.txt'%(snp_folder))
depth_cutoff = 6
min_good_alignment_samples = .4 # % of samples with bad alignment
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
    lineage_SNP.setdefault(lineage, [0, dict()]) # total SNPs on genes, gene SNPs
    for lines in open(genemut_file,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\t')
            if lines_set[6] != 'None':
                # only SNPs on a gene
                genename = lines_set[6]
                lineage_SNP[lineage][0] += 1
                lineage_SNP[lineage][1].setdefault(genename, [0,set(),[],0,0])# No. SNPs + genome set with SNPs, No. isolates, No. N SNPs, No. truncation SNPs
                lineage_SNP[lineage][1][genename][0] += 1
                lineage_SNP[lineage][1][genename][1].add(lines_set[5])
                alleles = lines_set[4]
                alleles_total = len(alleles)
                lineage_SNP[lineage][1][genename][2].append(
                    (len(alleles) - len([x for x in alleles if x == '-'])) / alleles_total)
                if lines_set[8] in ['N','NN']:
                    # N SNPs
                    lineage_SNP[lineage][1][genename][-1] += 1
                    if '*' in lines_set[9]:
                        # truncation
                        lineage_SNP[lineage][1][genename][-2] += 1
    return lineage_SNP

def load_clonal(clonal_file):
    clonal = dict()
    for lines in open(clonal_file,'r'):
        if not lines.startswith('cluster'):
            lines_set = lines.split('\n')[0].split('\t')
            lineage = lines_set[0]
            nonORFlength = int(lines_set[1])-int(lines_set[2])
            ORFlength = int(lines_set[2])
            clonal.setdefault(lineage,[nonORFlength,ORFlength])
    return clonal

def load_genes(assembly,allSNPscov):
    total_gene_length = 0
    gene_length = dict()
    depth_median = np.median(allSNPscov['avg_depth'])
    gene_copy_output = []
    for record in SeqIO.parse(assembly, 'fasta'):
        record_id = str(record.id)
        CHR = '_'.join(record_id.split('_')[:-1])
        record_des = str(record.description)
        startPOS = int(record_des.split(' # ')[1])-1000
        endPOS = int(record_des.split(' # ')[2])+1000
        allSNPscovsub = allSNPscov.loc[(allSNPscov['CHR'] == CHR) & (allSNPscov['POS'] <= endPOS) & (allSNPscov['POS'] >= startPOS), :]
        if allSNPscovsub.shape[0] > 0:
            record_seq_len = len(str(record.seq))
            gene_copy_number = np.median(allSNPscovsub['avg_depth'])/depth_median
            if gene_copy_number < max_depth_fold:
                gene_length.setdefault(record_id, [record_seq_len,gene_copy_number]) #gene length, depth fold = copy number
                total_gene_length += record_seq_len
                gene_copy_output.append('%s\t%s\n'%(record_id,gene_copy_number))
    f1 = open(assembly + '.copynum.txt','w')
    f1.write('gene_name\tcopy_number\n')
    f1.write(''.join(gene_copy_output))
    f1.close()
    return [gene_length,total_gene_length]

def find_clonal(lineage_short):
    if lineage_short in clonal:
        return clonal[lineage_short]
    else:
        return clonal.get(lineage_short.split('_')[0],[0,0])

def find_assemlby(lineage):
    return glob.glob('%s/%s.noHM.fasta'%(assembly_folder,lineage))

def pvalue_mutgene(SNP, ORFlength,gene_length,gene_set):
    SNP_genome_set_all = set()
    pvalueset = set()
    mut_rate = float(SNP)/ORFlength
    allgenepvalue = ['gene\tpvalue\n']
    for gene in gene_length:
        if gene in gene_set:
            SNP_gene,SNP_genome_set,num_strains_with_mut,truncation_gene,N_SNP_gene = gene_set.get(gene,[0,set(),[1],0,0])
            num_strains_with_mut = mean(num_strains_with_mut)
            this_gene_length,gene_copy_num = gene_length[gene]
            # NOT USED considering the prevalence of this gene among strains
            #pvalue = 1-poisson.cdf(SNP_gene - 1,mut_rate*this_gene_length*num_strains_with_mut)# greater than and equal to
            # treat all genes prevalence as 100%
            pvalue = 1 - poisson.cdf(SNP_gene - 1, mut_rate * this_gene_length)
            #k, n, p = [SNP_gene, this_gene_length, mut_rate * num_strains_with_mut]
            #pvalue = 1 - binom.cdf(k-1, n, p)   # greater than and equal to
            allsum_details.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.5f\t%s\n'%(lineage,gene,pvalue,SNP_gene,N_SNP_gene,truncation_gene,this_gene_length,gene_copy_num,num_strains_with_mut,mut_rate*1000,len(SNP_genome_set)))
            allgenepvalue.append('%s\t%s\n'%(gene,pvalue))
            if SNP_gene > 1 and len(SNP_genome_set) > 1 and pvalue <= pvalue_max:
                # potential PE
                pvalueset.add(pvalue)
                SNP_genome_set_all.add(len(SNP_genome_set))
        else:
            allgenepvalue.append('%s\t%s\n' % (gene, 1))
    pvalueset = list(pvalueset)
    pvalueset.sort(reverse=True)
    SNP_genome_set_all = list(SNP_genome_set_all)
    SNP_genome_set_all.sort()
    if len(pvalueset) > 0:
        foutput = open('%s/allgenes.poisson.pvalue.within.%s.details.txt' % (summary_folder,lineage), 'w')
        foutput.write(''.join(allgenepvalue))
        foutput.close()
    return [pvalueset,SNP_genome_set_all]

def mutation_sim(SNP, ORFlength,gene_length,pvalueset,SNP_genome_set_all):
    gene_num = []
    gene_length_set = []
    num_strains_with_mut = 1 - min_good_alignment_samples# lower bound
    allgenepvalue = ['simulation\tpvalue\n']
    mut_rate = float(SNP)/ORFlength
    for gene,this_gene_length in gene_length.items():
        this_gene_length,gene_copy_num = this_gene_length
        gene_num.append(gene)
        gene_length_set.append(this_gene_length)
    for i in range(0,simulation_round):
        gene_mut = random.choices(gene_num, weights=gene_length_set, k=SNP)
        allgenes_mut = set(gene_mut)
        allsim_gene = dict()
        for geneID in gene_num:
            if geneID in allgenes_mut:
                SNP_gene = gene_mut.count(geneID)
                this_gene_length,gene_copy_num = gene_length[geneID]
                #pvalue = 1 - poisson.cdf(SNP_gene - 1, mut_rate * this_gene_length * num_strains_with_mut) # greater than and equal to
                # treat all genes prevalence as 100%
                pvalue = 1 - poisson.cdf(SNP_gene - 1, mut_rate * this_gene_length)
                #k, n, p = [SNP_gene, this_gene_length, mut_rate * num_strains_with_mut]
                #pvalue = 1 - binom.cdf(k-1, n, p)  # greater than and equal to
                allgenepvalue.append('%s\t%s\n' % (i, pvalue))
                if SNP_gene > 1:
                    # FP PE
                    for SNP_gene_sub in range(2,SNP_gene):
                        allsim_gene.setdefault(SNP_gene_sub,[])
                        allsim_gene[SNP_gene_sub].append(pvalue)
            else:
                allgenepvalue.append('%s\t%s\n' % (i, 1))
        for SNP_gene_sub in SNP_genome_set_all:
            for pvalue in pvalueset:
                num_genes_pass_pvalue = len([x for x in allsim_gene.get(SNP_gene_sub,[]) if x <= pvalue])
                allsum.append('%s\t%s\t%s\t%s\t%s\n'%(lineage,SNP_gene_sub,pvalue,i,num_genes_pass_pvalue))
    if len(pvalueset) > 0:
        foutput = open('%s/allgenes.poisson.pvalue.within.simulation.%s.details.txt' % (summary_folder, lineage), 'w')
        foutput.write(''.join(allgenepvalue))
        foutput.close()
    return allsum

def depth_screening(allSNPscov):
    allSNPscov = allSNPscov[allSNPscov['POS'] > 100]  # not considering the depth of the first POS of each contig
    allSNPscov.index = range(0, allSNPscov.shape[0])
    # avg depth of non zero
    total_genomes = allSNPscov.shape[1] - 3
    for i in allSNPscov.index:
        non_zero_depth = [x for x in allSNPscov.iloc[i, 3:] if x > 0]
        good_depth = [x for x in allSNPscov.iloc[i, 3:] if x > depth_cutoff]
        if len(good_depth) >= min_good_alignment_samples * total_genomes:
            # good coverage, >= min_fraction_of_good_coverage isolates with >= min_depth depth
            allSNPscov.loc[i, 'avg_depth'] = np.mean(non_zero_depth)
        else:
            allSNPscov.loc[i, 'avg_depth'] = 0
    allSNPscov = allSNPscov[allSNPscov['avg_depth'] >= depth_cutoff] # good depth
    return allSNPscov

################################################### Main ########################################################
# load SNPs and genes with mutations
lineage_SNP = dict()
for genemut_file in allgenemut_file:
    print('processing %s'%(genemut_file))
    lineage_SNP = load_genemut(genemut_file,lineage_SNP)

# load non ORF length
clonal = load_clonal(clonal_file)

# simulation
allsum_details = []
allsum_details.append('lineage\tgene_name\tpvalue\tSNP_gene\tN_gene\ttruncation_gene\tgene_length\tcopy_number\tprevalence\tmut_rate_1kbp\tgenome_set\n')
allsum = []
allsum.append('lineage\tSNP_cutoff\tpvalue_cutoff\tsimulation_round\tFP\n')
for lineage in lineage_SNP:
    allSNPscov = pd.read_csv(os.path.join(MW_folder, '%s.cov.MW.txt') % (
        lineage), sep='\t')
    allSNPscov = depth_screening(allSNPscov)
    nonORFlength,ORFlength = find_clonal(lineage) # callable genome length
    SNP, gene_set = lineage_SNP[lineage]
    # load gene length
    assembly = find_assemlby(lineage)
    if assembly == []:
        print('no assembly for %s in %s' % (lineage, assembly_folder))
    else:
        assembly = assembly[0]
        # load gene length for all genes
        gene_length, total_gene_length = load_genes(assembly.replace('.fasta', '.fna'), allSNPscov)
        print('process %s %s SNPs %s all ORF length with coverage %s all ORF length with coverage and 1 copy' % (
        lineage, SNP, ORFlength, total_gene_length))
        ORFlength = total_gene_length  # not considering genes with >= 2 copies
        if SNP > 1 and gene_set != dict():  # SNPs on gene > 1
            lineage_new = lineage.replace('_clustercluster', '_CL').replace('_PB_', '_PaDi_').split('.donor')[0]
            # compute pvalue for all genes
            pvalueset, SNP_genome_set_all = pvalue_mutgene(SNP, ORFlength, gene_length, gene_set)
            print('finish compute poisson pvalue %s PE pvalue set %s SNP genome set %s' % (
            lineage, pvalueset, SNP_genome_set_all))
            if len(pvalueset) > 0:
                # simulation
                mutation_sim(SNP, ORFlength, gene_length, pvalueset, SNP_genome_set_all)
            print('finish simulation %s' % (lineage))

foutput = open('%s/allgenes.poisson.pvalue.within.txt'%(summary_folder), 'w')
foutput.write(''.join(allsum_details))
foutput.close()
foutput = open('%s/allgenes.poisson.pvalue.within.simulation.txt'%(summary_folder), 'w')
foutput.write(''.join(allsum))
foutput.close()
################################################### END ########################################################
