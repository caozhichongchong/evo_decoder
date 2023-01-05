################################################### END ########################################################
################################################### SET PATH ########################################################
# Filter results of WGS
import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
import statistics
import itertools
import random
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="a folder to store all output",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round3/',
                      metavar='output/')
required.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/donor_species/evo_decoder/',
                      metavar='scripts/')


################################################## Definition ########################################################
args = parser.parse_args()
# set up path
input_script = args.s + '/vcf_filter.sh'
output_dir = args.i + '/vcf/'
vcf_name = '.raw.vcf'
ref_filename = '.noHM.fasta'
Cov_dis_overall = 1000
MW_cov_dir = '%s/MV_cov'%(output_dir)

try:
    os.mkdir(MW_cov_dir)
except IOError:
    pass

f1 = open(os.path.join(input_script), 'w')
cmds = '#!/bin/bash\nsource ~/.bashrc\npy39\n'
f1.write(cmds)
f1.close()

def filter_vcfs(allvcfs,lineage):
    CHRPOS_withSNPs = os.path.join(output_dir, '%s.CHRPOS' % (lineage))
    try:
        f1 = open(CHRPOS_withSNPs,'r')
    except IOError:
        for vcf in allvcfs:
            snp_file = vcf.replace('.raw.vcf', '.flt.snp.vcf')
            os.system('cat %s | cut -f 1,2 >> %s' % (snp_file, CHRPOS_withSNPs))
        os.system('cat %s | sort | uniq > %s.unique'%(CHRPOS_withSNPs,CHRPOS_withSNPs))
    cmds = ''
    for vcf in allvcfs:
        try:
            f1 = open('%s.filter'%(vcf),'r')
        except IOError:
            cmds += 'bcftools filter -T %s.unique %s > %s.filter\n' % (CHRPOS_withSNPs, vcf, vcf)
    f1 = open(os.path.join(input_script), 'a')
    f1.write(cmds)
    f1.close()

def outputcovwindow(allvcfs,lineage):
    try:
        f1 = open(MW_cov_dir + '/%s.cov.MW.txt' % (lineage),'r')
    except IOError:
        cov_genome = dict()
        i = 0
        sample_len = len(allvcfs)
        for vcf in allvcfs:
            for lines in open(vcf, 'r'):
                if not lines.startswith('#'):
                    lines_set = lines.split('\n')[0].split('\t')
                    CHR = lines_set[0]
                    POS = int(lines_set[1])
                    if POS % Cov_dis_overall == 500:
                        CHRPOS = '%s\t%s' % (CHR, POS)
                        cov_genome.setdefault(CHRPOS, [0] * sample_len)
                        # subset this CHR
                        Subdepth_all = lines_set[9]
                        Subdepth = Subdepth_all.split(':')[-1].replace('\n', '').split(',')
                        total_sub_depth = sum([int(Subdepth_sub) for Subdepth_sub in Subdepth])
                        cov_genome[CHRPOS][i] = total_sub_depth
            i += 1
        cov_output = []
        allCHRPOS = list(cov_genome.keys())
        allCHRPOS.sort()
        for CHRPOS in allCHRPOS:
            temp_cov = cov_genome[CHRPOS]
            cov_output.append('%s\t%s\t%s\n' % (CHRPOS, statistics.mean(temp_cov),
                                                                '\t'.join(str(cov) for cov in temp_cov)))
        vcf_file_filtered = open(MW_cov_dir + '/%s.cov.MW.txt' % (lineage), 'w')
        vcf_file_filtered.write('CHR\tPOS\tavg_depth\t%s\n' % ('\t'.join([os.path.basename(vcf).split('.genome.')[1].split(vcf_name)[0] for vcf in allvcfs])) + ''.join(cov_output))
        vcf_file_filtered.close()

################################################### Main ########################################################
# run vcf filtering
all_vcf_file=glob.glob(os.path.join(output_dir,'*%s'%(vcf_name)))

# cluster donor species
all_vcf_file_set = dict()
for vcf_file in all_vcf_file:
    lineage = os.path.basename(vcf_file).split('.genome')[0]
    all_vcf_file_set.setdefault(lineage, [])
    all_vcf_file_set[lineage].append(vcf_file)

for lineage in all_vcf_file_set:
    allvcfs = all_vcf_file_set[lineage]
    allvcfs.sort()
    print('process',lineage)
    filter_vcfs(allvcfs, lineage)
    outputcovwindow(allvcfs,lineage)