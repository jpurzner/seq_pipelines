#!/srv/gs1/software/python/python-2.7/bin/python
import subprocess
import os
import re
from collections import defaultdict
import bradner_pipeline
import glob
import itertools 
import sys

# pipeline for processing ChIP-seq data
# for a given pool all fastq files should be stored in the same folder 

# map the fastq files


#=========================================================================
#=============================GLOBAL PARAMETERS===========================       
#=========================================================================  

fastqDelimiter = ','


#==========================================================================
#=======================FUNCTIONS==========================================


def load_meta(meta_file):
    '''
    loads a meta data table into a dictionary 
    the first column should be the desired filename 
    and the second the index for the fastq files 
    '''

    meta = {}    
    with open(meta_file) as mf:
        for line in mf:
            meta_line = line.split()
            meta[meta_line[1]] = meta_line[0]
            
    return meta
    
def get_fastq_names(meta):
    '''
    input: meta dictionary from load_meta
         key = index
         entry = filename (name of the bam files, etc)

    output: dictionary of fastq files for each index 
         key = index
         entry = fastq_files

    gets all fastq files with the index in the name and 
    stores them in a dictionary with the index as a key
    
    will display index files and warning if number of files
    different per index
    
    if the same filename is used for multiple indexes then 
    the fastq files associated with the index will be treated 
    as biological replicates and merged 
    
    '''
    print meta
    fastq_files = defaultdict(list)
    for f in os.listdir(os.curdir):
        # only check fastq files
        # 3/10/16 for some reason I had fastq. which prevented the reading 
        # of naked fastq files, changed
        if 'fastq' in f:
            # check each index                                            
            for index in meta.keys():
                if index in f:
                    fastq_files[meta[index]].append(f)
            
    # reverse meta to easily reference 
    metai = {v: k for k, v in meta.items()}
    # check the number of files per index 
    fastq_num = []                
    for name in fastq_files.keys():
        print 'index: ' + metai[name]
        print 'name: ' + name 
        print 'fastq files:' + '\n' + '\n'.join(fastq_files[name])
        curr_num = len(fastq_files[name])        
        if fastq_num and (curr_num not in fastq_num):
            print 'warning: number of fastq files per index not the same as others in directory'
        fastq_num.append(curr_num)
    
    return fastq_files


def genome_dict(machine = "scg3",  genome = 'mm9'): 
    '''
    this function holds the paths for genome and indexes for each of our systems
    '''

    if machine == 'scg3' and genome == 'mm9':
        genome_file = {
           'bowtie2': '/srv/gsfs0/projects/cho/bowtie2_indexes',
           'chr_size': '/srv/gsfs0/projects/cho/annotation/mm9_chr_size_noran.txt',
           'genome_chr': '/srv/gsfs0/projects/cho/genomes/mm9/chr',
           'mappability': '/srv/gsfs0/projects/cho/genomes/mm9/mappability/globalmap_k20tok54'    
        }
    elif machine == 'cho_oro' and genome == 'mm9':
        genome_file = {
            'bowtie2': '/tank/genomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index',
            'genome_chr': '/tank/genomes/Mus_musculus/UCSC/mm9/Sequence/Chromosomes',
            'genome': '/tank/genomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta',
            'mappability': '/tank/genomes/align2raw/umap/globalmap_k20tok54'
        }
    else: 
        print "error: machine not specified or no information given"
        print "please add paths and filenames to genome_dict"
    
    return genome_file


def bash_header(machine = 'scg3'):
    if machine == 'scg3':
        cmd = '#$ -l h_vmem=10G\n'
        cmd += '#$ -l h_rt=24:00:00\n'
        cmd += '#$ -w e\n'
        cmd += '#$ -cwd\n'
        cmd += '#$ -V\n'
        cmd += '#$ -pe shm 2\n'
    elif machine == 'cho_oro':
        cmd = ''

    return cmd 

def fastq_read_length():
    '''
    determines the size of the reads to optimize read trimming
    NEEDS TO BE COMPLETED, READ off compressed file
    PASS INFORMATION TO TRIMMOMATIC
    '''

    seq_list = []
    with open(fastq) as f:
        for line in itertools.islice(f, 2, 1002, 4):
            seq_list.append(line)

    print seq_list



def fastq_trim(trimmomaticString, fileNameDict, pairedEnd, sequence_length = 36, adapters = 'illumina_list.fa' ):
    '''
    creates the trimmomatic command                                                    
    '''
    trimlogFile = fileNameDict['trimlogFile']
    cmd = "touch %s\n" % (trimlogFile)
    adapters_file = os.path.dirname(trimmomaticString) + '/adapters/' + adapters
    
    min_length = sequence_length - 8 
    
    if pairedEnd:

        tempFastqFile1 = fileNameDict['tempFastqFile_1']
        tempFastqFile2 = fileNameDict['tempFastqFile_2']
        
    else:
        tempFastqFile = fileNameDict['tempFastqFile']
        tempFastqFileTrim = fileNameDict['tempFastqFileTrim']
        
        if isinstance(tempFastqFile, list): # fastqFile is a list          
            for i, fastq in enumerate(tempFastqFile):
                cmd += "java -Xmx6g -Xms3g -Djava.awt.headless=true -jar %s SE -phred33 %s %s  ILLUMINACLIP:%s:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:%s >> %s 2>&1\n" % (trimmomaticString, fastq, tempFastqFileTrim[i], adapters_file,  str(min_length), trimlogFile)
        else:
            cmd += "java -Xmx6g -Xms3g -Djava.awt.headless=true -jar %s SE -phred33 %s %s  ILLUMINACLIP:%s:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:%s >> %s 2>&1\n" % (trimmomaticString, fastq, tempFastqFileTrim, adapters_file,  str(min_length), trimlogFile)
                        
    return cmd


def rmdupBamCmd(samtoolsString,fileNameDict):

    '''
    uses smatools to sort the bam
    '''

    logFile = fileNameDict['logFile']
    tempSortedBamFile = fileNameDict['tempSortedBamFile']    
    tempRmdupBamFile = fileNameDict['tempRmdupBamFile']
    cmd = "echo \"Remove duplicate log\" >> %s\n " % (logFile)   
    cmd += "%s rmdup -s  '%s' '%s' >> %s 2>&1" % (samtoolsString,tempSortedBamFile,tempRmdupBamFile, logFile)
    return cmd


def bam2tagAlign(samtoolsString, fileNameDict):
    
    '''
    generate tagAlign file for align2rawsignal
    reads modified to fit the unique mappability tracks 
    excuse the confusing awk code
    '''

    tempRmdupBamFile = fileNameDict['tempRmdupBamFile'] 
    temptagAlignFile = fileNameDict['temptagAlignFile']
    
    # innumerable escape characters to make allow this string to proceed without errors
    cmd  = "%s view -F 0x0204 %s | gawk 'BEGIN{OFS=\"\\t\"}{if (and($2,16) > 0) {print $3,($4-1),($4-1+length($10)),\"N\",\"1000\",\"-\"} else {print $3,($4-1),($4-1+length($10)),\"N\",\"1000\",\"+\"} }' | gawk 'BEGIN {OFS=\"\\t\"} {if ((($3-$2)>20) && (($3-$2)<54)) print $0; else if ((($3-$2)>54) && ($6==\"+\")) print $1, $2, $2+53, $4, $5, $6; else if ((($3-$2) > 54) && ($6==\"-\")) print $1, $3-53, $3, $4, $5, $6}' > %s" % (samtoolsString, tempRmdupBamFile, temptagAlignFile)

    return cmd


def wiggler(matlabString, wigglersrc,  fileNameDict, genome_files):
    
    '''
    generate matlab commands to produce normalized tracks using 
    align2rawsignal https://code.google.com/p/align2rawsignal/
    written by Anshul Kundaje
    '''

    temptagAlignFile = fileNameDict['temptagAlignFile']
    wigglerFile = fileNameDict['wigglerFile']
    wigglerLog = fileNameDict['wigglerLog']
    logFile = fileNameDict['logFile']

    path_command = 'addpath(\'%s;\')' % (wigglersrc)
    mat_command = "align2rawsignal('-i=%s', '-s=%s', '-u=%s', '-o=%s', '-of=bg', '-mm=12', '-l=200', '-v=%s');" % (temptagAlignFile, genome_files['genome_chr'], genome_files['mappability'], wigglerFile, wigglerLog)
    
    cmd = 'touch wiggler.m\n'
    cmd += 'echo \"%s\" >> wiggler.m\n' % (path_command)
    cmd += 'echo \"%s\"  >> wiggler.m\n' % (mat_command)
    cmd += 'echo "quit;" >> wiggler.m\n'
    cmd += '%s -nojvm -nodisplay -nosplash -nodesktop -r wiggler.m\n' % (matlabString)
    cmd += 'rm wiggler.m\n'
    #cmd += 'echo \"Align2rawsignal log:\" >> %s\n' % (logFile)
    #cmd += 'cat %s >> %s\n' % (wigglerLog, logFile)
    return cmd

def bedgraph2tdf(tdfString, fileNameDict, genome = 'mm9'):
    ''' 
    converts bedgraph files to tdf files for visualization with IGV
    '''

    wigglerFile = fileNameDict['wigglerFile']
    tdfFile = fileNameDict['tdfFile']

    cmd = "java -Xmx6g -Xms500m  -Djava.awt.headless=true -jar %s toTDF -z 5  %s %s %s" % (tdfString, wigglerFile, tdfFile, genome)  

    return cmd

def rmtagAlingCmd(fileNameDict):
    '''    
    remove the tagAlign file
    '''

    temptagAlignFile = fileNameDict['temptagAlignFile']
    cmd = "/bin/rm -f '%s'" % (temptagAlignFile)
    return cmd

def rmbamFile(fileNameDict):
    tempBamFile = fileNameDict['tempBamFile']
    cmd = "/bin/rm -f '%s'" % (tempBamFile)
    return cmd


def map_qual(samtoolsString, mapqString, fileNameDict):
    '''
    generate the histogram of mapping qualities
    using MAPQ_hist.R 
    '''

    tempBamFile = fileNameDict['tempBamFile']
    logFile = fileNameDict['logFile']
    mapqFile = fileNameDict['mapqFile']

    cmd = "%s view %s | cut -f 5 | sort -n | uniq -c > %s\n" % (samtoolsString, tempBamFile, mapqFile)
    cmd += "Rscript %s %s\n" % (mapqString, mapqFile)
    cmd += "echo \"Mapping quality summary\" >> %s\n" % (logFile)
    cmd += "cat %s >> %s\n" % (mapqFile, logFile)
    cmd += "rm %s\n" % (mapqFile)
    return cmd

def bedgraph_check(fileNameDict): 
    '''
    count the number of chromosomes in each bedgraph file
    '''
    wigglerFile = fileNameDict['wigglerFile']
    bedchrFile = fileNameDict['bedchrFile']    
    logFile = fileNameDict['logFile']
    cmd = "echo \"#chromosomes that are found in the bedgraph file\" > %s\n" % (bedchrFile)
    cmd += "echo \"#if chromosomes are missing then the wiggler didn't complete\"  > %s\n" % (bedchrFile)
    cmd += "cat %s | cut -f 1 | uniq > %s\n" % (wigglerFile, bedchrFile)
    cmd += "echo \"chromosomes in align2rawsignal output\" >>  %s\n" % (logFile)
    cmd += "cat %s | tr \'\n\' \' \' >> %s\n" % (bedchrFile, logFile)
    cmd += "echo \"\n\" >> %s \n" % (logFile)
    cmd += "echo \"number of chromosomes in bedgraph\" >> %s \n" % (logFile)
    cmd += "wc -l %s >> %s\n" % (bedchrFile, logFile)
    return cmd

def gz_bedgraph(fileNameDict):
    wigglerFile = fileNameDict['wigglerFile']
    cmd = "gzip %s\n" % (wigglerFile)
    return cmd 


def bash_maker(fileNameDict, outputFolder, uniqueID, pairedEnd, mismatchN, genome_files, exe_dict, genome = 'mm9', trim = True, machine = 'scg3', readlen = 36):
    
    #open the bashfile to write to                          
    bashFileName = "%s%s_bwt.sh" % (outputFolder,uniqueID)
    bashFile = open(bashFileName,'w')
    
    #shebang                  
    bashFile.write('#!/bin/bash\n')

    #add machine specific header 
    cmd = bash_header()
    bashFile.write(cmd+'\n')

    #make temp directory      
    cmd = 'mkdir %s' % (fileNameDict['tempFolder'])
    bashFile.write(cmd+'\n')
    
    #make final directory
    cmd = 'mkdir %s' % (fileNameDict['finalFolder'])
    bashFile.write(cmd+'\n')
    
    #extract fastq    
    cmd = bradner_pipeline.extractFastqCmd(fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')
    
    # create the not_complete file
    cmd = 'touch %s' % (fileNameDict['not_complete'])
    bashFile.write(cmd+'\n')

    #call fastqc                                                     
    cmd = bradner_pipeline.runFastQC(exe_dict['fastqc'],fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')
    
    # trim reads
    if trim:
        cmd = fastq_trim(exe_dict['trimmomatic'], fileNameDict, pairedEnd)
        bashFile.write(cmd+'\n')

        #rerun fastqc 
        cmd = bradner_pipeline.runFastQC(exe_dict['fastqc'],fileNameDict,pairedEnd,trim)
        bashFile.write(cmd+'\n')

    #call bowtie    
    cmd = bradner_pipeline.bowtieCmd(exe_dict['bowtie2'],mismatchN,genome_files['bowtie2'],fileNameDict, genome = genome, pairedEnd = pairedEnd, trim = trim)
    bashFile.write(cmd+'\n')

    #remove temp fastq                                 
    cmd = bradner_pipeline.removeTempFastqCmd(fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #remove temp trimmed fastq 
    if trim:
        cmd = bradner_pipeline.removeTempFastqCmd(fileNameDict,pairedEnd, trim = trim)
        bashFile.write(cmd+'\n')

    #generate a bam                                   
    cmd = bradner_pipeline.generateTempBamCmd(exe_dict['samtools'],fileNameDict)
    bashFile.write(cmd+'\n')

    #check mapq
    cmd = map_qual(exe_dict['samtools'], exe_dict['MAPQ_hist'], fileNameDict)
    bashFile.write(cmd+'\n')
    
    #change into the temp directory
    cmd = bradner_pipeline.changeTempDir(fileNameDict)
    bashFile.write(cmd+'\n')

    #sort the bam
    cmd = bradner_pipeline.sortBamCmd(exe_dict['samtools'],fileNameDict)
    bashFile.write(cmd+'\n')

    #index    
    cmd = bradner_pipeline.indexBamCmd(exe_dict['samtools'],fileNameDict)
    bashFile.write(cmd+'\n')

    #remove duplicates
    cmd = rmdupBamCmd(exe_dict['samtools'], fileNameDict) 
    bashFile.write(cmd+'\n')

    #index the removed duplicates bam 
    cmd = bradner_pipeline.indexBamCmd(exe_dict['samtools'],fileNameDict, rmdup = True)
    bashFile.write(cmd+'\n')

    #generate wiggler file
    cmd = bam2tagAlign(exe_dict['samtools'], fileNameDict)
    bashFile.write(cmd+'\n')

    cmd = wiggler(exe_dict['matlab'], exe_dict['wiggler_src'], fileNameDict, genome_files)
    bashFile.write(cmd+'\n')

    #generate tdf file
    cmd = bedgraph2tdf(exe_dict['igvtools'], fileNameDict, genome = genome)
    bashFile.write(cmd+'\n')
    
    #check wiggler by counting the chromosomes in the bedgraph file
    cmd = bedgraph_check(fileNameDict)
    bashFile.write(cmd+'\n')

    #compress the bedgraph file  
    cmd = gz_bedgraph(fileNameDict)
    bashFile.write(cmd+'\n')

    #remove sam       
    cmd = bradner_pipeline.rmSamCmd(fileNameDict)
    bashFile.write(cmd+'\n')
    
    #remove tagAlign files 
    cmd = rmtagAlingCmd(fileNameDict)
    bashFile.write(cmd+'\n')
    
    #remove unsorted bam file
    cmd = rmbamFile(fileNameDict)
    bashFile.write(cmd+'\n')

    #mv bams                                        
    cmd = bradner_pipeline.mvBamCmd(fileNameDict)
    bashFile.write(cmd+'\n')

    #cleanup        
    cmd = bradner_pipeline.rmTempFiles(fileNameDict)
    bashFile.write(cmd+'\n')

    # create the not_complete file                                               
    cmd= 'rm %s' % (fileNameDict['not_complete'])
    bashFile.write(cmd+'\n')

    bashFile.close()

    return bashFileName 

def run_script(bashFileName, machine):
    if machine == 'scg3':
        cmd = "qsub %s" % (bashFileName)
        subprocess.call(cmd, shell = True)
    elif machine == 'cho_oro':
        cmd = "nohup bash %s &" % (bashFileName) 
        subprocess.call(cmd, shell = True)
    else:
        sys.exit( "error: incorrect machine specified")


def main():
    '''
    main 
    '''

    import argparse
    import random
    
    # obtain command line arguments                                            
    parser = argparse.ArgumentParser(description='batch pipeline for mapping ChIP-seq data')
    parser.add_argument('-m','--meta', help='table containing filename<tab>index', required=True)
    parser.add_argument('-o','--output', help='directory to create new files in, default is current working directory', default='',  required=False)
    parser.add_argument('-d','--directory', help='directory storing fastq files, default current working directory', default='',  required=False)
    parser.add_argument('-c','--machine', help='speficies the paths for genome.fa, chr sizes and bowtie2 index for a given system', default='scg3', choices=['scg3','cho_oro'], required=False)
    parser.add_argument('-g','--genome', help='specifies the genome to be used for mapping', default='mm9', required=False)
    parser.add_argument('-p','--paired', help='flag specifies paired end, without flag assumes single end reads', action='store_true', required=False)
    parser.add_argument('-l','--readlen', help='specify the length of reads, default 36 bp', default=36,  required=False)
    parser.add_argument('--no_trim', help='flag to prevent trimming reads using Trimmomatic prior to mapping, default trimming', action='store_true', required=False)
    parser.add_argument('-t','--temp', help='specifies the temporary folder path, otherwise will use current working directory', default='', required=False)
    parser.add_argument('--bowtie_mismatch', help='number of mismatches for bowtie2', default=0, required=False)
    parser.add_argument('-n','--norun', action='store_true', required=False, help='only generate the bash files but do not run them')

    
    args = parser.parse_args()

    
    #======================================================================#
    #============== PROCESS INPUT ARGUMENTS ===============================#

    # process meta data and validate existence of index 
    meta = load_meta(args.meta)
    fastq_dict = get_fastq_names(meta)
    machine = args.machine
    genome = args.genome
    genome_files = genome_dict(machine = machine, genome = genome) 
    pairedEnd = args.paired
    
    trim = not args.no_trim
    

    readlen = args.readlen
    mismatchN = args.bowtie_mismatch
    norun = args.norun 

    # setup temporary directory
    # if no temp directory use current working directory
    if args.temp == '':
        tempParentFolder = os.getcwd() + '/'
    else:
        tempParentFolder = args.temp
        if not tempParentFolder.endswith('/'):
            tempParentFolder = tempParentFolder + '/'

    # get current directory
    if args.directory == '':
        fastq_directory = os.getcwd() + '/'
    else:
        fastq_directory = args.directory
        if not fastq_directory.endswith('/'):
            fastq_directory = fastq_directory + '/'

    # setup output directory
    # if no output directory use current working directory              
    if args.output == '':
        outputFolder = os.getcwd() + '/'
    else:
        outputFolder = args.ouput
        if not outputFolder.endswith('/'):
            outputFolder = outputFolder + '/'
    
    print "test test"
    #=====================================================================#     
    #=================== GET PATH FOR EXECUTABLES ========================#

    # find executables
    exe_list = ['bowtie2', 'samtools', 'fastqc', 'matlab']
    exe_dict = {}
    for e in exe_list:
        which_str = 'which' + ' ' + e 
        
        res = subprocess.Popen(['which', e], stdout=subprocess.PIPE)
        output, err = res.communicate()
        exe_dict[e] = output.strip('\n')
        
    # find programs that use JRE and R scripts
    search_path = os.environ.get('PATH')
    
    path_list = search_path.split(':')
    for path_name in path_list:
        trimmomatic_temp = glob.glob(path_name + '/trimmomatic*')
        igvtools_temp = glob.glob(path_name + '/igvtools.jar')
        MAPQ_hist_temp = glob.glob(path_name + '/MAPQ_hist.R')

        if trimmomatic_temp:
            trimmomatic = ''.join(trimmomatic_temp)
               
        if igvtools_temp:
            igvtools = ''.join(igvtools_temp)
                    
        if MAPQ_hist_temp:
            MAPQ_hist = ''.join(MAPQ_hist_temp)
            
    if trimmomatic:
        exe_dict['trimmomatic'] = trimmomatic
    else:
        print 'error: timmomatic not autodetected in path'
        sys.exit('add to path, install, or ommit trimming')
    if igvtools:
        exe_dict['igvtools'] = igvtools    
    else: 
        print 'error: igvtools not autodetected in path'
        sys.exit('add to path, install, or omit tdf creation')

    if MAPQ_hist:
        exe_dict['MAPQ_hist'] = MAPQ_hist
    else:
        print 'error: MAPQ_hist.R not autodetected in path'
        sys.exit('add to path, install, or omit MAPQ plotting')

    exe_dict['wiggler_src'] = '/srv/gsfs0/projects/cho/programs/align2rawsignal/src'
    #print exe_dict

    #======================================================================#
    #=============== GENERATE BASH SCRIPTS & EXECUTE=======================#

    print fastq_dict
    for name in fastq_dict.keys():
        fastqFile = fastq_dict[name] 
        uniqueID = name

        # generate temporary string for folder and file names  
        tempString = '_%s' % str(random.randint(1,10000))
        fileNameDict = bradner_pipeline.makeFileNameDict(fastqFile, fastq_directory, tempString,tempParentFolder,outputFolder, uniqueID = uniqueID, pairedEnd = pairedEnd)
        
        bashFileName = bash_maker(fileNameDict, outputFolder, uniqueID, pairedEnd, mismatchN, genome_files, exe_dict, genome = genome, trim = trim, machine = machine, readlen = readlen)
        
        if norun: 
            print "bash scripts generated without running"
        else: 
            run_script(bashFileName, machine)

        




    

if __name__ == "__main__":
    main()
