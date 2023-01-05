import glob
import os
import statistics
# set up path
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument("-p",
                      help="path file of all fastqs",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/allgenome_path.txt',
                      metavar='allgenome_path.txt')
required.add_argument("-o",
                      help="output folder",
                      type=str, default='/scratch/users/anniz44/genomes/donor_species/vcf_round3/',
                      metavar='output/')
required.add_argument("-s",
                      help="a folder to store all scripts",
                      type=str, default='/scratch/users/anniz44/scripts/1MG/donor_species/evo_decoder/',
                      metavar='scripts/')
################################################## Definition ########################################################
args = parser.parse_args()
outputfolder = args.o
input_script = args.s
input_script_vcf = '%s/mapping_WGS'%(input_script)
working_dir = os.getcwd()
try:
    os.mkdir(outputfolder)
except IOError:
    pass
try:
    os.mkdir(outputfolder + '/co-assembly')
except IOError:
    pass
try:
    os.mkdir(outputfolder + '/vcf')
except IOError:
    pass
try:
    os.mkdir(outputfolder + '/co-assembly-species')
except IOError:
    pass
try:
    os.mkdir(outputfolder + '/vcf-species')
except IOError:
    pass
try:
    os.mkdir(input_script)
except IOError:
    pass
try:
    os.mkdir(input_script_vcf)
except IOError:
    pass

genomesize_range = [0.9,1.1]  #gemove genomes out of genome size * X-Y

def run_vcf_WGS(files,files2,database,tempbamoutput):
    # generate code
    cmds = ''
    try:
        f1 = open('%s.sorted.bam' % (tempbamoutput),'r')
    except IOError:
        cmds += 'bowtie2' + ' --threads %s -x %s -1 %s -2 %s |%s view -@ %s -S -b >%s.bam\n%s sort -@ %s %s.bam -o %s.sorted.bam\n%s index -@ %s %s.sorted.bam\n' % (
            40, database, files, files2, 'samtools', 40,
            tempbamoutput, 'samtools', 40, tempbamoutput, tempbamoutput, 'samtools', 40,
            tempbamoutput)
        cmds += 'rm -r %s.bam %s.bam.bai\n' % (tempbamoutput, tempbamoutput)
    return [cmds, '%s.sorted.bam' % (tempbamoutput)]

def merge_sample(database,vcfoutput,allsam):
    cmds = ''
    try:
        f1 = open('%s.raw.vcf' % (vcfoutput), 'r')
    except FileNotFoundError:
        cmds += '%s mpileup --threads %s -a FMT/ADF,FMT/ADR,FMT/AD -q30 -Ou -B -d3000 -f %s %s | %s call -c -Ov --threads %s > %s.raw.vcf\n' % (
            'bcftools', 40, database,
            ' '.join(allsam), 'bcftools', 40, vcfoutput)
        cmds += 'rm -r %s\n'%(' '.join(allsam))
    try:
        f1 = open('%s.flt.snp.vcf' % (vcfoutput))
    except FileNotFoundError:
        cmds += '%s view -H -v snps %s.raw.vcf > %s.flt.snp.vcf \n' % (
            'bcftools', vcfoutput, vcfoutput)
    return cmds

def check_genome_size(fastas,genomename,fastq):
    allgenomesize = []
    for fasta in fastas:
        try:
            filesize = int(os.path.getsize(fasta))
            allgenomesize.append(filesize)
        except FileNotFoundError:
            allgenomesize.append(0)
    # filter out genomesize that's above the size cutoff
    avg_genome_size = statistics.median(allgenomesize)
    min_genome_size = genomesize_range[0]*avg_genome_size
    max_genome_size = genomesize_range[1]*avg_genome_size
    qualify_fastqs = [fastq[i] for i in range(0, len(allgenomesize)) if
                     allgenomesize[i] <= max_genome_size and allgenomesize[i] >= min_genome_size]
    qualify_genomes = [genomename[i] for i in range(0, len(allgenomesize)) if
                     allgenomesize[i] <= max_genome_size and allgenomesize[i] >= min_genome_size]
    not_qualified = [x for x in genomename if x not in qualify_genomes]
    print(len(qualify_fastqs),len(not_qualified))
    print(avg_genome_size,[allgenomesize[genomename.index(x)]/avg_genome_size for x in not_qualified])
    return [qualify_fastqs,qualify_genomes]

def load_path(pathfile):
    alllineage = dict()
    newresult = dict()
    for lines in open(pathfile,'r'):
        if not lines.startswith('genome'):
            lines_set = lines.split('\n')[0].split('\t')
            genome,lineage,fastq,fasta = lines_set
            if fastq != '' and fasta != '':
                alllineage.setdefault(lineage,[[],[],[]])
                alllineage[lineage][0].append(genome)
                alllineage[lineage][1].append(fastq)
                alllineage[lineage][2].append(fasta)
                newresult.setdefault(genome,lines)
    return [alllineage,newresult]

def subset(file1,file2,output_name1,output_name2):
    cmds = 'head -400000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'tail -400000 %s >> %s\n' % (
        file1, output_name1)
    cmds += 'head -400000 %s >> %s\n' % (
        file2, output_name2)
    cmds += 'tail -400000 %s >> %s\n' % (
        file2, output_name2)
    return cmds

def runspades(file1,file2,output_name):
    temp_output = output_name.split('.fasta')[0]
    cmds = '%s --careful -1 %s -2 %s -o %s --threads 40 --memory 100 --cov-cutoff 7\n' % \
            ('spades.py',file1, file2, temp_output)
    cmds += 'mv %s/scaffolds.fasta %s\n' % (temp_output, output_name)
    cmds += 'rm -r %s %s %s\n' % (temp_output, file1, file2)
    return cmds

def runprodigal(fasta):
    try:
        f1 = open('%s.fna'%(fasta.split('.fasta')[0]),'r')
        return ''
    except IOError:
        return 'prodigal -q -i %s -d %s.fna -a %s.faa\n'%(fasta,fasta.split('.fasta')[0],fasta.split('.fasta')[0])

def runprokka(output_folder,fasta):
    cmds = ''
    try:
        f1 = open('%s.gff'%(fasta),'r')
    except IOError:
        cmds += 'prokka --kingdom Bacteria --force --outdir %s %s\n'%(output_folder, fasta)
        cmds += 'mv %s/*.gff %s.gff\n'%(output_folder,fasta)
        cmds += 'rm -r %s\n'%(output_folder)
    return cmds

alllineage,newresult = load_path(args.p)
# process each lineage
qualified = ['genome\tlineage\tfastq\tfasta\n']
allspecies_count = dict()
for lineage in alllineage:
    print('processing %s'%(lineage))
    species = lineage.split('_')[0]
    genomename,fastq,fasta = alllineage[lineage]
    # check genome size
    fastq,genomename = check_genome_size(fasta,genomename,fastq)
    allspecies_count.setdefault(species,[])
    if len(genomename) >= 3:
        # run mapping
        cmds = '#!/bin/bash\nsource ~/.bashrc\npy39\n'
        cmds_assembly = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        assembly_output = outputfolder + '/co-assembly/%s.fasta'%(lineage)
        assembly_output_noHM = outputfolder + '/co-assembly/%s.noHM.fasta' % (lineage)
        allspecies_count[species].append(assembly_output_noHM.replace('.fasta','.fna'))
        merge_fastq1 = outputfolder + '/co-assembly/%s_1.fastq'%(lineage)
        merge_fastq2 = outputfolder + '/co-assembly/%s_2.fastq' % (lineage)
        qualified += [newresult[x] for x in genomename]
        cmds += 'bowtie2-build %s %s\n' % (assembly_output_noHM, assembly_output_noHM)
        for i in range(0,len(genomename)):
            genome = genomename[i]
            fastq_file = fastq[i]
            fastq_file2 = fastq_file.replace('_1.fastq', '_2.fastq')
            cmds_assembly += subset(fastq_file, fastq_file2, merge_fastq1, merge_fastq2)
            try:
                f1 = open('%s.raw.vcf' % (outputfolder + '/vcf/%s.genome.%s' % (lineage,genome)), 'r')
            except IOError:
                results = run_vcf_WGS(fastq_file, fastq_file2,
                                      assembly_output_noHM,
                                      outputfolder + '/vcf/%s.genome.%s' % (lineage,genome))
                cmds += results[0]
                cmds += merge_sample(assembly_output_noHM, outputfolder + '/vcf/%s.genome.%s' % (lineage,genome), [results[1]])
        cmds_assembly += runspades(merge_fastq1, merge_fastq2,
                                   assembly_output)
        cmds_prokka = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        cmds_prokka += runprokka(assembly_output + '_prokka',assembly_output_noHM.replace('.fasta','.fna'))
        cmds_prodigal = '#!/bin/bash\nsource ~/.bashrc\npy37\n'
        cmds_prodigal += runprodigal(assembly_output_noHM)
        f1 = open(os.path.join(input_script_vcf, '%s.vcf.sh' % (lineage)), 'w')
        try:
            checkassembly = open(assembly_output,'r')
        except IOError:
            f1.write(cmds_assembly)
        try:
            checkassembly = open(assembly_output_noHM, 'r')
        except IOError:
            cmds_HM = '#!/bin/bash\nsource ~/.bashrc\nmodule add c3ddb/glibc/2.14\nblastn -subject %s -query %s -outfmt 6 -out %s.blasn.txt -max_target_seqs 1000 -num_threads 40 -perc_identity %s \n'%(
                assembly_output,assembly_output,assembly_output,90)
            cmds_HM += 'py37\npython %s/remove_homologous.py -i %s -blastn %s.blasn.txt\n'%(working_dir,assembly_output,assembly_output)
            f1.write(cmds_HM)
        f1.write(cmds_prodigal)
        f1.write(cmds_prokka)
        f1.write(cmds)
        f1.close()

# process each species

for species in allspecies_count:
    allassembly = allspecies_count[species]
    if len(allassembly) > 1:
        cmds_HM = '#!/bin/bash\nsource ~/.bashrc\nmodule add c3ddb/glibc/2.14\n'
        # at least 2 lineages
        print('processing %s'%(species))
        assembly1 = allassembly[0]
        for assembly2 in allassembly:
            cmds_HM += 'blastn -subject %s -query %s -outfmt 6 -out %s_%s.blasn.txt -max_target_seqs 1000 -num_threads 40 -perc_identity %s -qcov_hsp_perc %s \n'%(
                assembly1,assembly2,assembly2.split('.noHM')[0],os.path.basename(assembly1).split('.noHM')[0],90,85)
        f1 = open(os.path.join(input_script_vcf, '%s.blast.sh' % (species)), 'w')
        f1.write(cmds_HM)
        f1.close()

# sum all codes
f1 = open(os.path.join(input_script, 'allWGS.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.vcf.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allWGS.sh'%(input_script))

f1 = open(os.path.join(input_script, 'allblast.sh'), 'w')
f1.write('#!/bin/bash\nsource ~/.bashrc\n')
for sub_scripts in glob.glob(os.path.join(input_script_vcf, '*.blast.sh')):
    sub_scripts_name = os.path.split(sub_scripts)[-1]
    f1.write('jobmit %s %s\n' % (sub_scripts, os.path.split(sub_scripts)[-1]))

f1.close()
print('please run: sh %s/allblast.sh'%(input_script))

f1 = open('/scratch/users/anniz44/genomes/donor_species/vcf_round2/merge/details/allgenome_path_qualified.txt','w')
f1.write(''.join(qualified))
f1.close()