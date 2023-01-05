# sum up aligned region ORF length and non-ORF length -> cluster_length_new.py
import os,glob
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

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
assemblyfolder = args.i + '/co-assembly/'
fasta = '.noHM.faa'
allcovfiles = glob.glob(
'%s/*.txt'%(args.i + '/vcf/MV_cov/'))
print(allcovfiles)
try:
    os.mkdir(args.i + '/mergevcf/')
except IOError:
    pass

outputfile = '%s/clonal_genelength_new.txt'%(args.i + '/mergevcf/')
max_fraction_ambigious_samples = .4 #If more than % of the samples have ambiguous NT, discard the candidate location
min_fraction_of_good_coverage = 1-max_fraction_ambigious_samples
Cov_dis_overall = 1000 # calculate coverage per 1000 bp
min_depth = 6 #Remove candidate locations have lower than this depth

def load_covfile(covfile,fnafile):
    allgenes = dict()
    for record in SeqIO.parse(fnafile, 'fasta'):
        record_id = str(record.id)
        chr = '_'.join(record_id.split('_')[:-1])
        allgenes.setdefault(chr,{})
        record_descriptionset = str(record.description).split(' # ')
        allgenes[chr].setdefault(int(record_descriptionset[1]),int(record_descriptionset[2])-int(record_descriptionset[1])+1) # start, gene length
    covresult = [0,0,0] # genomic POS covered, ORF-region POS, number of genes
    for lines in open(covfile,'r'):
        if not lines.startswith('CHR'):
            lines_set = lines.split('\n')[0].split('\t')
            if len([int(x)>=min_depth for x in lines_set[3:]]) >= min_fraction_of_good_coverage*(len(lines_set)-3):
                # >= min_fraction_of_good_coverage isolates with >= min_depth depth
                # good coverage
                CHR, POS = lines_set[:2]
                POS = int(POS)
                if POS > 1:
                    covresult[0] += Cov_dis_overall # genomic POS covered
                    # check genes
                    allgeneschr = allgenes.get(CHR,{})
                    checked_genes = []
                    for genestart in allgeneschr:
                        if genestart<=POS and genestart >= POS - Cov_dis_overall:
                            # genes within this moving window region
                            covresult[1] += allgeneschr[genestart] # ORF-region POS covered
                            covresult[2] += 1 # number of genes covered
                            checked_genes.append(genestart)
                        elif genestart > POS:
                            break
                    for genestart in checked_genes:
                        # remove genes that already checked
                        allgeneschr.pop(genestart)
                    allgenes[chr] = allgeneschr
    return [str(x) for x in covresult]

alloutput = ['cluster\tgenome_coverage\tORF_region_coverage\tnum_genes_covered\n']
for covfile in allcovfiles:
    covfilename = os.path.split(covfile)[-1]
    lineagename = covfilename.split('.cov.MW.txt')[0]
    fnafile = glob.glob('%s/%s%s'%(assemblyfolder,
                                   lineagename,
                                   fasta))
    if len(fnafile) > 0:
        print('process %s'%(covfile))
        covresult = load_covfile(covfile,fnafile[0])
        alloutput.append('%s\t%s\n'%(lineagename,
                         '\t'.join(covresult)))
        print(covresult)

f1 = open(outputfile,'w')
f1.write(''.join(alloutput))
f1.close()

