#!/usr/bin/python
#given a fastq location, create a .sh script to run mapping through generation of sorted bam

'''
The MIT License (MIT)

Copyright (c) 2013 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#===================================================================
#========================MODULES AND DEPENDENCIES===================
#===================================================================

import string
import random
import os
import sys 

#===================================================================
#==========================GLOBAL PARAMETERS========================
#===================================================================

#command arguments
bowtieString = 'bowtie2'
samtoolsString = 'samtools'
fastqcString = '/usr/local/FastQC/fastqc'
fastqDelimiter = ','

#tempParentFolder = '/mnt/d0-0/share/bradnerlab/projects/anna/BOWTIE_TEMP/'

#===================================================================
#=============================FUNCTIONS=============================
#===================================================================


def stripExtension(fileName):

    '''
    tries to strip the extension of a filename
    can strip .tar.gz, .txt, .fastq,.gz,.zip
    modified to handle lists of filenames 
    '''
    extensionList = ['.tar.gz','.tar', '.txt','.fastq','.fasta','.gz','.zip', '.fastq.gz']

    if isinstance(fileName, list):
        stripName = []
        for f in fileName: 
           for extension in extensionList:
               f = f.replace(extension,'')
           stripName.append(f)
    else: 
        for extension in extensionList:
            stripName  = fileName.replace(extension,'')
 
    return stripName
    

def makeFileNameDict(fastqFile, fastqDir,tempString,tempParentFolder,finalParentFolder,linkFolder = '', mismatchN = 0, uniqueID='',pairedEnd = False):

    '''
    creates a dictionary w/ all filenames

    modified to take lists of fastq files and to make names for fastq trimming
    '''

    fileNameDict = {}
    
    fastqName = uniqueID
    fileNameDict['fastqName'] = uniqueID

    #make the temp Folder and final folder 
    tempFolder = tempParentFolder + 'bwt_' + fastqName + tempString + '/'
    fileNameDict['tempFolder'] = tempFolder

    finalFolder = finalParentFolder + fastqName + '/'


    if pairedEnd:
        # check that there are even number of entries 
        # the order of the fastq files maintains the index values  
        print(len(fastqFile))
        if (len(fastqFile) % 2 == 0):  
            if (len(fastqFile) == 2):    
                fileNameDict['fastqFile_1'] = fastqDir + fastqFile[0]
                fileNameDict['fastqFile_2'] = fastqDir + fastqFile[1]
                
                tempFastqFile1 = tempFolder + fastqName + '_1.rawFastq'
                fileNameDict['tempFastqFile_1'] = tempFastqFile1
                tempFastqFileTrim = tempFolder + fastqName + '_1_trim' + '.rawFastq'
                fileNameDict['tempFastqFileTrim_1'] = tempFastqFileTrim
               
                tempFastqFile2 = tempFolder + fastqName + '_2.rawFastq'
                fileNameDict['tempFastqFile_2'] = tempFastqFile2
                tempFastqFileTrim = tempFolder + fastqName + '_2_trim' + '.rawFastq'
                fileNameDict['tempFastqFileTrim_2'] = tempFastqFileTrim

            else: 
                fileNameDict['fastqFile_1'] = []
                fileNameDict['fastqFile_2'] = []
                fileNameDict['tempFastqFile_1'] = []
                fileNameDict['tempFastqFile_2'] = []
                fileNameDict['tempFastqFileTrim_1'] = []
                fileNameDict['tempFastqFileTrim_2'] = []

                # loop through the list of fastq files and add to appropriate group
                for i in range(0, len(fastqFile), 2):
                    fileNameDict['fastqFile_1'].append(fastqDir + fastqFile[i])
                    fileNameDict['fastqFile_2'].append(fastqDir + fastqFile[i+1])

                    fileNameDict['tempFastqFile_1'].append(tempFolder 
                          + stripExtension(fastqFile[i]) + '.rawFastq')
                    fileNameDict['tempFastqFile_2'].append(tempFolder 
                          + stripExtension(fastqFile[i+1]) + '.rawFastq')

                    fileNameDict['tempFastqFileTrim_1'].append(tempFolder
                          + stripExtension(fastqFile[i]) + '_trim' +  '.rawFastq')
                    fileNameDict['tempFastqFileTrim_2'].append(tempFolder 
                          + stripExtension(fastqFile[i+1]) + '_trim' + '.rawFastq')
        else: 
            sys.exit("error: paired end called but odd number of fastq files given")
  
    # for both single end and paired end files
    # aslo assign names to the same fastqFile used by single end to ease 
    # fastqc and removal of temp files 
    fileNameDict['fastqFile'] = [fastqDir + s  for s in fastqFile]       
    tempFastqFile = [tempFolder + stripExtension(f) + '.rawFastq' for f in fastqFile]
    tempFastqFileTrim = [tempFolder + stripExtension(f) + '_trim'  + '.rawFastq' for f in fastqFile]
    fileNameDict['tempFastqFile'] = tempFastqFile
    fileNameDict['tempFastqFileTrim'] = tempFastqFileTrim

    print(fileNameDict)


        # modified to handle list of fastq files
        #if isinstance(fastqFile, list): # fastqFile is a list
        #    tempFastqFile = [tempFolder + stripExtension(f) + '.rawFastq' for f in fastqFile]
        #    tempFastqFileTrim = [tempFolder + stripExtension(f) + '_trim'  + '.rawFastq' for f in fastqFile]
        #    fileNameDict['tempFastqFile'] = tempFastqFile
        #    fileNameDict['tempFastqFileTrim'] = tempFastqFileTrim
        #else:
        #    tempFastqFile = tempFolder + fastqName + '.rawFastq'
        #    tempFastqFileTrim = tempFolder + fastqName + '_trim' + '.rawFastq'
        #    fileNameDict['tempFastqFile'] = tempFastqFile
        #    fileNameDict['tempFastqFileTrim'] = tempFastqFileTrim
    
    tophatBam = tempFolder + 'tophat/accepted_hits.bam'
    fileNameDict['tophatBam'] = tophatBam

    logFile = tempFolder + fastqName + '.log'
    fileNameDict['logFile'] = logFile

    mapqFile = tempFolder + fastqName + '_mapq.txt'
    fileNameDict['mapqFile'] = mapqFile

    trimlogFile = tempFolder + fastqName + '.trim_log'
    fileNameDict['trimlogFile'] = trimlogFile
    
    tempSamFile = tempFolder + fastqName + '.sam'
    fileNameDict['tempSamFile'] = tempSamFile

    tempBamFile = tempFolder + fastqName + '.bam'
    fileNameDict['tempBamFile'] = tempBamFile

    finalBamFile = finalFolder + fastqName + '.bam'
    fileNameDict['finalBamFile'] = finalBamFile

    tempSortedBamFile = finalFolder + fastqName + '_sort.bam'
    fileNameDict['tempSortedBamFile'] = tempSortedBamFile

    tempRmdupBamFile = finalFolder + fastqName + '_sort_rmdup.bam'
    fileNameDict['tempRmdupBamFile'] = tempRmdupBamFile

    temptagAlignFile = finalFolder + fastqName + '.tagAlign'
    fileNameDict['temptagAlignFile'] = temptagAlignFile

    wigglerFile = finalFolder + fastqName + '_norm.bedgraph'
    fileNameDict['wigglerFile'] = wigglerFile
    
    wigglerLog = finalFolder + fastqName + '_norm.log'
    fileNameDict['wigglerLog'] = wigglerLog

    bedchrFile = finalFolder + fastqName + '_bedgraph_chr.txt'
    fileNameDict['bedchrFile'] = bedchrFile

    tdfFile = finalFolder + fastqName + '_norm.tdf'
    fileNameDict['tdfFile'] = tdfFile

    groupHeader = tempFolder + fastqName 
    fileNameDict['groupHeader'] = groupHeader

    fileNameDict['finalFolder'] = finalFolder

    fileNameDict['linkFolder'] = linkFolder
    
    not_complete = finalFolder + fastqName + '.not_complete'
    fileNameDict['not_complete'] = not_complete 

    return fileNameDict


def extractFastqCmd(fileNameDict,pairedEnd = False):

    '''
    creates a command to extract/copy the fastq to a temp location
    '''
    if pairedEnd:
        fastqList = []
        fastqFile1 = fileNameDict['fastqFile_1']
        tempFastqFile1 = fileNameDict['tempFastqFile_1']
        fastqFile2 = fileNameDict['fastqFile_2']
        tempFastqFile2 = fileNameDict['tempFastqFile_2']

        if isinstance(fastqFile1, list): # fastqFile is a list 
            fastqList1 = [list(i) for i in zip(fastqFile1, tempFastqFile1)]
            fastqList2 = [list(i) for i in zip(fastqFile2, tempFastqFile2)]
            fastqList = fastqList1 + fastqList2
        else: 
            fastqList.append(fastqFile1)
            fastqList.append(tempFastqFile1)
            fastqList.append(fastqFile2)
            fastqList.append(tempFastqFile2)
    else:
        fastqFile = fileNameDict['fastqFile']
        tempFastqFile = fileNameDict['tempFastqFile']
        fastqList = [list(i) for i in zip(fastqFile, tempFastqFile)]

    print(fastqList)
    cmdList = []

    for [fastqFile,tempFastqFile] in fastqList:  # initially used to handle paired end 
    #there are 3 possibilities, a gzipped, tarballed, or naked fastq
        if string.lower(fastqFile).count('tar.gz') == 1:
            cmd = "tar --strip-components 5 --to-stdout -xzvf %s > %s" % (fastqFile,tempFastqFile)
        elif string.lower(fastqFile).count('tar') == 1:
            cmd = "tar -xzvf %s > %s" % (fastqFile,tempFastqFile)
        elif string.lower(fastqFile.split('.')[-1]) == 'gz':
            cmd = 'cp %s %s.gz\n' % (fastqFile,tempFastqFile)
            cmd+= 'gunzip %s.gz' % (tempFastqFile)
        elif string.lower(fastqFile.split('.')[-1]) == 'fastq.gz':
            cmd = 'cp %s %s.fastq.gz\n' % (fastqFile,tempFastqFile)
            cmd+= 'gunzip %s.fastq.gz' % (tempFastqFile)
        else:
            cmd = 'cp %s %s' % (fastqFile,tempFastqFile)
    
        cmdList.append(cmd)

    fullCmd = string.join(cmdList,'\n')
    #print fullCmd

    return fullCmd

def runFastQC(fastqcString,fileNameDict,pairedEnd = False, trim = False):

    '''
    cmd to run fastqc
    modified to accept lists of fastq files and run on trimmed data 
    

    if pairedEnd: ##needs to be fixed 
        fastqName = fileNameDict['fastqName']
        tempFastqFile1 = fileNameDict['tempFastqFile_1']
        tempFastqFile2 = fileNameDict['tempFastqFile_2']
        finalFolder = fileNameDict['finalFolder']

        if finalFolder[-1] != '/':
            finalFolder+='/'
        finalFolder1 = finalFolder + '%s_1_fastqc' % (fastqName)
        finalFolder2 = finalFolder + '%s_2_fastqc' % (fastqName)
        
        cmd = 'mkdir %s\n' % (finalFolder1)
        cmd += 'mkdir %s\n' % (finalFolder2)

        if isinstance(tempFastqFile, list): # fastqFile is a list
            for f in tempFastqFile:
                cmd += '%s -o %s %s\n' % (fastqcString,finalFolder1,tempFastqFile1)
                cmd += '%s -o %s %s\n' % (fastqcString,finalFolder1,tempFastqFile1)
        else: 
            cmd += '%s -o %s %s\n' % (fastqcString,finalFolder1,tempFastqFile1)
            cmd += '%s -o %s %s' % (fastqcString,finalFolder2,tempFastqFile2)
    '''   
    fastqName = fileNameDict['fastqName']
    finalFolder = fileNameDict['finalFolder']
    if finalFolder[-1] != '/':
        finalFolder+='/'

    # if read trimming selected then re-run on trimmed fastq file
    if trim: 
        tempFastqFile = fileNameDict['tempFastqFileTrim']
        finalFolder += '%s_fastqc_trim' % (fastqName)
    else:
        tempFastqFile = fileNameDict['tempFastqFile']    
        finalFolder += '%s_fastqc' % (fastqName)

    cmd = 'mkdir %s\n' % (finalFolder)
    if isinstance(tempFastqFile, list): # fastqFile is a list
        for f in tempFastqFile:
            cmd += '%s -o %s %s\n' % (fastqcString,finalFolder,f)
    else:
        cmd += '%s -o %s %s' % (fastqcString,finalFolder,tempFastqFile)
    return cmd

def bowtieCmd(bowtieString,mismatchN,bowtieIndex,fileNameDict, genome = 'mm9', pairedEnd=False, trim = True):

    '''
    creates the bowtie command call
        -x maximum fragment length 
        -N mismatch 
        -x genome (bowtie indexes) 
    '''
    
    logFile = fileNameDict['logFile']    
    tempSamFile = fileNameDict['tempSamFile']
    cmd = "export BOWTIE2_INDEXES=%s\n" % (bowtieIndex)
    
#calling bowtie
    if pairedEnd:
        if trim: 
            tempFastqFile1 = fileNameDict['tempFastqFileTrim_1']
            tempFastqFile2 = fileNameDict['tempFastqFileTrim_2']
        else: 
            tempFastqFile1 = fileNameDict['tempFastqFile_1']
            tempFastqFile2 = fileNameDict['tempFastqFile_2']
        if isinstance(tempFastqFile1, list): # fastqFile is a list
            comma_fastq1 = ','.join(tempFastqFile1)    
            comma_fastq2 = ','.join(tempFastqFile2)    
            cmd += "%s -p 4 -X 2000 -N %s -x %s -1 %s -2 %s -S %s >> %s 2>&1" % (bowtieString,mismatchN, genome,comma_fastq1,comma_fastq2,tempSamFile, logFile)
        else:
            cmd += "%s -p 4 -X 2000 -N %s -x %s -1 %s -2 %s -S %s >> %s 2>&1" % (bowtieString,mismatchN, genome,tempFastqFile1,tempFastqFile2,tempSamFile, logFile)

    else:
        if trim:
            bowtieFastq = fileNameDict['tempFastqFileTrim']
        else:
            bowtieFastq = fileNameDict['tempFastqFile']


        if isinstance(bowtieFastq, list): # fastqFile is a list
            comma_fastq = ','.join(bowtieFastq)
            cmd += "%s -p 4 -N %s -x %s -U %s -S %s >> %s 2>&1" % (bowtieString,mismatchN,genome,comma_fastq,tempSamFile, logFile)
        else:
            cmd += "%s -p 4 -N %s -x %s -U %s -S %s >> %s 2>&1" % (bowtieString,mismatchN,genome,bowtieFastq,tempSamFile, logFile)
    return cmd


def removeTempFastqCmd(fileNameDict,pairedEnd = False, trim = False):

    '''
    removes the temp fastq
    '''

    if trim: 
        tempFastqFile = fileNameDict['tempFastqFileTrim']
    else:
        tempFastqFile = fileNameDict['tempFastqFile']
    if isinstance(tempFastqFile, list): # fastqFile is a list
        cmd = ''
        for f in tempFastqFile:
            cmd += '/bin/rm -f %s\n' % (f)     
    else:
        cmd = '/bin/rm -f %s' % (tempFastqFile)
    
    return cmd


def generateTempBamCmd(samtoolsString,fileNameDict):

    '''
    uses samtools to convert the sam to a bam
    '''

    tempSamFile = fileNameDict['tempSamFile']
    tempBamFile = fileNameDict['tempBamFile']
    cmd = "%s view -bS '%s' > '%s'" % (samtoolsString,tempSamFile,tempBamFile)
    
    return cmd

#change into temp directory

def changeTempDir(fileNameDict):

    '''
    changes into the temp directory
    '''

    tempFolder = fileNameDict['tempFolder']
    cmd = "cd %s" % (tempFolder)
    return cmd

#sort
def sortBamCmd(samtoolsString,fileNameDict):

    '''
    uses smatools to sort the bam
    '''
    finalBamFile = fileNameDict['finalBamFile']
    tempSortedBamFile = os.path.splitext(fileNameDict['tempSortedBamFile'])[0]
    
    cmd = "%s sort '%s' '%s'" % (samtoolsString,finalBamFile,tempSortedBamFile)
    return cmd


#index
def indexBamCmd(samtoolsString,fileNameDict, rmdup = False):

    '''
    uses samtools to index the bam
    '''
    
    if rmdup:
        BamFile = fileNameDict['tempRmdupBamFile']
    else:
        BamFile = fileNameDict['tempSortedBamFile']
        
    cmd = "%s index '%s'" % (samtoolsString, BamFile)

    return cmd

def rmSamCmd(fileNameDict):

    '''
    remove the sam
    '''

    tempSamFile = fileNameDict['tempSamFile']
    cmd = "/bin/rm -f '%s'" % (tempSamFile)
    return cmd

def mvBamCmd(fileNameDict):

    '''
    moves and renames the bam w/o the temp string
    '''

    groupHeader = fileNameDict['groupHeader']
    finalFolder = fileNameDict['finalFolder']

    cmd = "mv %s* %s" % (groupHeader,finalFolder)
    
    return cmd

def linkBamCmd(fileNameDict):

    '''
    moves and renames the bam w/o the temp string
    '''

    groupHeader = fileNameDict['groupHeader']
    finalFolder = fileNameDict['finalFolder']
    linkFolder = fileNameDict['linkFolder']

    cmd = "ln %s%s* %s" % (finalFolder,groupHeader,linkFolder)
    
    return cmd


def rmTempFiles(fileNameDict):

    '''
    removes the temp folder
    '''
    tempFolder = fileNameDict['tempFolder']
    
    cmd = "/bin/rm -r '%s'" % (tempFolder)
    return cmd

#===================================================================
#=============================MAIN==================================
#===================================================================

def main():

    '''
    main run function
    '''

    from optparse import OptionParser

    usage = "usage: %prog [options] -f [FASTQFILE] -g [GENOME] -u [UNIQUEID] -o [OUTPUTFOLDER]"
    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-f","--fastq", dest="fastq",nargs = 1, default=None,
                      help = "Enter the full path of a fastq file to be mapped")
    parser.add_option("-g","--genome",dest="genome",nargs =1, default = None,
                      help = "specify a genome, options are hg18 or mm9 right now")
    parser.add_option("-u","--unique",dest="unique",nargs =1, default = None,
                      help = "specify a uniqueID")
    parser.add_option("-o","--output",dest="output",nargs =1, default = None,
                      help = "Specify an output folder")


    #optional arguments
    parser.add_option("-N","--mismatch",dest="mismatchN",nargs =1, default = 0,
                      help = "Specify 0 or 1 for allowed mismatches")
    parser.add_option("--link-folder",dest="linkFolder",nargs =1, default = None,
                      help = "Specify a folder to symlink the bam")
    parser.add_option("-p","--paired",dest="paired",action='store_true',default = False,
                      help = "Flag for paired end data")



    (options,args) = parser.parse_args()

    if not options.fastq or not options.genome or not options.unique or not options.output:
        parser.print_help()
        exit()


    #retrive the arguments
    fastqFile = options.fastq
    genome = string.lower(options.genome)
    uniqueID = options.unique
    outputFolder = options.output

    #retrieve optional arguments
    mismatchN = options.mismatchN
    if options.linkFolder:

        linkFolder = options.linkFolder
    else:
        linkFolder =''
    pairedEnd = options.paired

    #get the bowtie index
    bowtieDict = {
        'mm9':'/raider/index/mm9/Bowtie2Index/genome',
        'hg19':'/raider/index/hg19/Bowtie2Index/genome',
        'hg18':'/grail/genomes/Homo_sapiens/human_gp_mar_06_no_random/bowtie/hg18'
        }

    bowtieIndex = bowtieDict[string.lower(genome)]

    #get the temp string
    tempString = '_%s' % str(random.randint(1,10000))
    
    fileNameDict = makeFileNameDict(fastqFile,genome,tempString,tempParentFolder,outputFolder,linkFolder,mismatchN,uniqueID,pairedEnd)

    #open the bashfile to write to
    bashFileName = "%s%s_bwt.sh" % (outputFolder,uniqueID)
    bashFile = open(bashFileName,'w')

    #shebang
    bashFile.write('#!/usr/bin/bash\n')

    #make temp directory
    cmd = 'mkdir %s' % (fileNameDict['tempFolder'])
    bashFile.write(cmd+'\n')

    #extract fastq
    cmd = extractFastqCmd(fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #call fastqc
    cmd =runFastQC(fastqcString,fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #call bowtie
    cmd = bowtieCmd(bowtieString,mismatchN,bowtieIndex,fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #remove temp fastq
    cmd = removeTempFastqCmd(fileNameDict,pairedEnd)
    bashFile.write(cmd+'\n')

    #generate a bam
    cmd = generateTempBamCmd(samtoolsString,fileNameDict)
    bashFile.write(cmd+'\n')

    #change into the temp directory
    cmd = changeTempDir(fileNameDict)
    bashFile.write(cmd+'\n')

    #sort the bam
    cmd = sortBamCmd(samtoolsString,fileNameDict)
    bashFile.write(cmd+'\n')

    #index
    cmd = indexBamCmd(samtoolsString,fileNameDict)
    bashFile.write(cmd+'\n')

    #remove sam
    cmd = rmSamCmd(fileNameDict)
    bashFile.write(cmd+'\n')

    #mv bams
    cmd = mvBamCmd(fileNameDict)
    bashFile.write(cmd+'\n')

    #link bams
    if options.linkFolder:
        cmd = linkBamCmd(fileNameDict)
        bashFile.write(cmd+'\n')

    #cleanup
    cmd = rmTempFiles(fileNameDict)
    bashFile.write(cmd+'\n')


    bashFile.close()

    print "Wrote mapping command to %s" % (bashFileName)
if __name__ == "__main__":
    main()

