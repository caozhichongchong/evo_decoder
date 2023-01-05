import glob
import os
import statistics
from Bio import SeqIO
from Bio.Seq import Seq
# set up path
import argparse
############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-i",
                      help="path of genome",
                      type=str, default='.',
                      metavar='input/genome.fasta')
optional.add_argument('-blastn',
                      help="blast result path",
                      metavar="genome.blastn.txt",
                      action='store', default='genome.blastn.txt', type=str)

############################################ Functions ##############################################
args = parser.parse_args()
# homologous cutoff
length_cutoff = 1000 # 1000 bp for homologous region
identity_cutoff = 90 # 90% identity for homologous region
CHR_length_cutoff = 1000 # minimum contig lengths for reference genome

def length_CHR(CHR):
    try:
        total_length = CHR.split('size')[1]
    except IndexError:
        total_length = CHR.split('length_')[1].split('_cov')[0]
    return int(total_length)

def remove_homologous(genome):
    print('removing homologous regions for genome %s' % (genome))
    HM_region = dict()
    for lines in open(args.blastn, 'r'):
        lines_set = lines.split('\n')[0].split('\t')
        CHR1, CHR2, identity, length_HM = lines_set[0:4]
        CHR2_start, CHR2_end = lines_set[8:10]
        if CHR1 != CHR2 and float(identity) >= identity_cutoff and float(length_HM) >= length_cutoff:
            # homologous region cutoff
            HM_region.setdefault(CHR2, [])
            HM_region[CHR2].append([int(CHR2_start) - 1, int(CHR2_end) - 1])
    print(HM_region)
    Newgenome = []
    for record in SeqIO.parse(genome, 'fasta'):
        record_id = str(record.id)
        record_seq = str(record.seq)
        if len(record_seq) >= CHR_length_cutoff:
            if record_id not in HM_region:
                Newgenome.append('>%s\n%s\n' % (record_id.split('_cov')[0], record_seq))
            else:
                CHRset = HM_region[record_id]
                record_seq = list(record_seq)
                for CHR_start, CHR_end in CHRset:
                    record_seq[CHR_start: (CHR_end + 1)] = 'N'*(CHR_end - CHR_start)
                record_seq = ''.join(record_seq)
                Newgenome.append('>%s\n%s\n' % (record_id.split('_cov')[0], record_seq))
    f1 = open('%s.noHM.fasta' % (genome.split('.fasta')[0]), 'w')
    f1.write(''.join(Newgenome))
    f1.close()

############################################ Main ##############################################
remove_homologous(args.i)
