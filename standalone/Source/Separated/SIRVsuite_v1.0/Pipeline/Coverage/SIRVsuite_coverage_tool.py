import numpy as np
import pandas as pd
import pysam as ps
from scipy.stats import norm
import re
import os
import glob 
import cairo
from Bio import SeqIO
import argparse
import matplotlib.colors as colors
import pyBigWig as bwig
import drawer

class SIRVsuite_coverage:
    
    def __init__(self, annotationFile, outputPath="/", geneList = [],
                UTRmode=False, transitionLengths=False, sampleNames = [], experimentName = "", UTRlength = 200,
                average=False, fasta=None, mix=None, uniqueMap=False, sort = False, strandSpecificLib = False, createBigWig = False):
        
        # Initialize default values
        
        self.UTRmode = UTRmode
        self.UTRlength = UTRlength
        self.annotationFile = annotationFile
        self.genes = geneList
        self.mixPath = mix
        self.strands = ["+","-"]
        self.err = ""
        self.outPath = outputPath
        self.showTransitionLengths = transitionLengths
        self.fastaPath = fasta
        self.sort = sort
        self.uniqueMap = uniqueMap
        self.average = average
        self.nonStrandedLib = strandSpecificLib
        self.createCovFile = False
        self.canProceed = True
        self.CODCalc = True
        self.createBigWig = createBigWig
        self.sampleNames = sampleNames
        self.experimentName = experimentName
        
        ## IF NO GENE SPECIFIED, TAKING ALL FROM THE ANNOTATION
        if ( len(geneList) == 0):
            self.allGenes = True
        else:
            self.allGenes = False
            
    def tagsExtract(self,featureNames,line):

        main_colNames = ["seqname","source","feature","start","end","score","strand","frame"]

        if(any(item in featureNames for item in main_colNames)):    
            pattern_main_columns = re.compile('\s*(.*?)\t')
            main_col_result = pattern_main_columns.findall(line)
            if (len(main_col_result)!=8):
                print ("It seems format of gtf file is wrong..")
                return

        matches = []
        for name in featureNames:

            if ( name in main_colNames ):
                    index = main_colNames.index(name)
                    matches.append(main_col_result[index])
            else:
                pattern = re.compile(name+'\s"(.*?)";')
                result = pattern.search(line)

                if ( result == None ):
                    matches.append("NaN")
                else:
                    matches.append(result.group(1))

        return matches
    
    def extractFromPath(self,filePath,feature):
        if ( feature == "ext" ):
            pattern = ".*\.(.*)"
        elif ( feature == "parentDir" ):
            pattern = ".*\/(.*)\/"
        elif ( feature == "fileName" ):
            pattern = "(?!.*\/)(.*)\."
        elif ( feature == "all" ):
            pattern = "(\/.*\/)(.*)\.(.*)"
        else:
            print ("selected feature is not defined")
            return

        result = re.compile(pattern)
        result = result.search(filePath)
        if ( result != None ):
            result = result.group(1)
            return result
        else:
            print ("There was a problem with path format")
            return

    def readGtf(self,fileName):
        with open(fileName, 'r') as file_object: 
            requiredTags = ["seqname","feature","start","end","strand","gene_id","transcript_id"]
            df = pd.DataFrame(columns=requiredTags)

            for line in file_object:
                tags = self.tagsExtract(requiredTags,line)
                if ( tags[1] != "exon" ):
                    continue
                df = df.append(pd.Series(tags,index=df.columns),ignore_index=True)
            try:  
                df["start"] = df["start"].astype(int)
                df["end"] = df["end"].astype(int)
            except ValueError:
                print ("GTF file in wrong format!")
                return

            return df
        
    def deleteTemporaryFile(self):
        os.system('rm -r '+ self.outPath+"temp")
            
    def reduceGTF(self):
        if ( os.path.isdir(self.outPath+"temp/") ):
            self.deleteTemporaryFile()
            
        os.makedirs(self.outPath+"temp/")
        
        command = ""
        
        for gene in self.genes:
            command += gene + "\|"
            
        command = command[0:-2]
        
        os.system('grep "'+command+'" '+self.annotationFile+' > '+self.outPath+'"temp/reduced.gtf"')
        
    def readAnotationFile(self,annotationFile,fileType):
        
        if ( fileType == "gtf" ):
            with open(annotationFile, 'r') as file_object: 
                
                requiredTags = ["seqname","feature","start","end","strand","gene_id","transcript_id"]
                extracted_DF = pd.DataFrame(columns=requiredTags)

                for line in file_object:
                    
                    comment_test = re.match("\s*#.*",line)
                    if comment_test is not None:
                        continue
                    
                    tags = self.tagsExtract(requiredTags,line)
                    
                    if (tags[1] != "exon"):
                        continue
                    extracted_DF = extracted_DF.append(pd.Series(tags,index=extracted_DF.columns),ignore_index=True)
                try:  
                    extracted_DF["start"] = extracted_DF["start"].astype(int)
                    extracted_DF["end"] = extracted_DF["end"].astype(int)
                except ValueError:
                    print ("GTF file in wrong format!")
            
        elif ( (fileType == "fasta") | (fileType == "fa") ):
            fasta_chrom_ID = np.empty(0)
            fasta_seq = np.empty(0)

            extracted_DF = pd.DataFrame(data={'chromosome':[],'sequence':[]})

            for record in SeqIO.parse(self.fastaPath, "fasta"):
                extracted_DF = extracted_DF.append(pd.Series([record.id,record.seq],index=extracted_DF.columns),ignore_index=True)
        
        elif ( fileType == "mix" ):
            
            mixes = pd.read_csv(self.mixPath,header=None,sep="\t",engine="python",comment="#")
            extracted_DF = pd.DataFrame(mixes)
            extracted_DF = extracted_DF.replace('"|;','',regex=True)
            extracted_DF[2] = extracted_DF[2].astype(float)

            if( np.shape(extracted_DF)[1] != 4 ):
                print ("Mix file isn't in correct format")
                return
        
        else:
            print ("uknown annotation type: "+fileType)
            extracted_DF = []

        return (extracted_DF)

    def annot2UTR(self,annotation_df):
        
        UTR_annotation = annotation_df[0:0] # Creating an empty data frame with the same structure as annotation data frame

        for strand in self.strands:
            transcripts = np.unique(annotation_df[annotation_df["strand"]  ==  strand]["transcript_id"])

            for transcript in transcripts:

                annot_transcript = annotation_df[annotation_df["transcript_id"]  ==  transcript]
                UTR_length = self.UTRlength

                if ( strand  ==  '+'):
                    annot_transcript = annot_transcript.sort_values(by=["start"], axis = 0, ascending = False).reset_index(drop=True)
                    exon_lengths = annot_transcript["end"] - annot_transcript["start"]

                    remainder = UTR_length
                    UTR_data_frame = annot_transcript

                    for index, length in enumerate(exon_lengths):
                        remainder = remainder - length
                        if (remainder < 0):
                            UTR_data_frame.loc[index,"start"] = UTR_data_frame.loc[index,"start"] - remainder
                            UTR_data_frame = UTR_data_frame.loc[0:index,:]
                            break

                elif ( strand  ==  '-'):
                    annot_transcript = annot_transcript.sort_values(by=["start"], axis = 0, ascending = True).reset_index(drop=True)
                    exon_lengths = annot_transcript["end"] - annot_transcript["start"]

                    remainder = UTR_length
                    UTR_data_frame = annot_transcript

                    for index, length in enumerate(exon_lengths):
                        remainder = remainder - length
                        if ( remainder < 0 ):
                            UTR_data_frame.loc[index,"end"] = UTR_data_frame.loc[index,"end"] + remainder
                            UTR_data_frame = UTR_data_frame.loc[0:index,:]
                            break
                
                UTR_annotation = UTR_annotation.append(UTR_data_frame)

        UTR_annotation = UTR_annotation.reset_index(drop=True)

        return UTR_annotation
        
    def process_annotation(self):
        ## This functions reads annotations files from input array and selects data for specific gene 
        ## and only exon information
        
        print ("Reading annotation files:")           
        
        ## LOADING FASTA
        if ( self.fastaPath != None ):
            if ( not os.path.isfile(self.fastaPath) ):
                print ("fasta file doesn't exist.. please try again")
                return
            else:
                print ("loading "+self.fastaPath)
                fileType = self.extractFromPath(self.fastaPath,"ext")
                self.fastaData = self.readAnotationFile(self.fastaPath,fileType)
        
        ## LOADING MIX
        if ( self.mixPath != None ):
            if ( not os.path.isfile(self.mixPath) ):
                print ("mix file doesn't exist.. please try again")
                return
            else:
                print ("loading "+self.mixPath)
                fileType = self.extractFromPath(self.mixPath,"ext")
                self.mix = self.readAnotationFile(self.mixPath,fileType)
        else:
            self.mix = pd.DataFrame(np.array([["test","test","test","test"]]))
            print ("Mixfile not provided.. transcripts for visualization taken from annotation..")
            
        ## LOADING GTF
        if ( not os.path.isfile(self.annotationFile) ):
            print ("annotation file doesn't exist.. please try again")
            return
        else:
            dataFrames = {"whole":{}}
            fileType = self.extractFromPath(self.annotationFile,"ext")
            print ("loading "+self.annotationFile)

            if ( self.UTRmode ):
                dataFrames["UTR_only"] = {}

            if ( not self.allGenes ):
                self.reduceGTF()
                dataFrame = self.readAnotationFile(self.outPath+"temp/reduced.gtf", fileType)
                self.genes = np.unique(dataFrame["gene_id"])
                self.deleteTemporaryFile()
            else:
                dataFrame = self.readAnotationFile(self.annotationFile,fileType)

            for index in np.arange(0,len(self.genes)):

                if ( self.annotationFile[-3:] == "gtf" ):

                    selectedDF = dataFrame[(dataFrame["gene_id"] == self.genes[index])]

                    ## LIMITING GAPDH TO THE MOST DOMINANT TRANSCRIPT ##
                    if ( self.genes[index]  ==  "ENSG00000111640" ):
                        selectedDF = dataFrame[(dataFrame["transcript_id"] == "ENST00000229239")]

                    if ( np.shape(selectedDF)[0] == 0 ):
                        print ("Gene: "+self.genes[index]+" not found in the annotation..")
                        return 
                    else:
                        dataFrames["whole"][self.genes[index]] = selectedDF
                        if ( self.UTRmode ):
                            dataFrames["UTR_only"][self.genes[index]] = self.annot2UTR(selectedDF)   

                elif ( self.annotationFile[0][-3:] == "bed" ):
                    print ("bed file yet to be developed")
        
        self.selectedAnnotations = dataFrames


    def calculate_real_coverage_BigWig(self, bigwig_path):
        
        realCov = {}
        print ("Calculating real coverage..")

        for file_index, bigwig_path in enumerate(bigwig_path):

            bigwig_path_full = {'+': bigwig_path + "_plus.bw",
                                '-': bigwig_path + "_minus.bw"}

            sampleName = self.sampleNames[file_index]

            realCov[sampleName] = {}
            realCov["genePositions"] = {}

            for gene in self.genes:
                realCov[sampleName][gene] = {}
                realCov[sampleName][gene][self.strands[0]]= {}
                realCov[sampleName][gene][self.strands[1]] = {}

                if (self.UTRmode):
                    gene_annotation = self.selectedAnnotations["whole"][gene]
                else:
                    gene_annotation = self.selectedAnnotations["whole"][gene]

                chromosomeNames = np.unique(gene_annotation["seqname"]).astype(str)
                start_genes = np.array(gene_annotation["start"])
                end_genes = np.array(gene_annotation["end"])

                if ( len(chromosomeNames) > 1 ):
                    print ("didn't expect multiple chromosome for one gene..")
                    return
                
                start_genes = np.array(gene_annotation["start"])
                end_genes = np.array(gene_annotation["end"])

                start = np.min(start_genes) - 1
                end = np.max(end_genes)
                
                for strand in self.strands:
                    file = bwig.open(bigwig_path_full[strand])
                    realCov[sampleName][gene][strand] = file.values(chromosomeNames[0],start,end)
            
        for strand in self.strands:
            for gene in self.genes:
                realCov['genePositions'][gene] = np.arange(start,end+1)

        if ( self.average ):
            realCov["avg"] = {}
            for gene in self.genes:
                realCov["avg"][gene]={}
                for strand in self.strands:
                    realCov["avg"][gene][strand]={}
                    realCov["avg"][gene][strand] = np.zeros(len(realCov["genePositions"][gene]))
                    for sampleName in self.sampleNames:
                        realCov["avg"][gene][strand] += realCov[sampleName][gene][strand]
                    realCov["avg"][gene][strand] /= len(self.sampleNames)

            self.sampleNames.append("avg")
                
        self.realCoverage = realCov 

    def calculate_real_coverage(self, bamFileList):
        # this method returns two numpy arrays of coverage as first output and reverse complemented coverage as second
        # specify the path of the bam file and array of reference genes
        # if you want to select all of reference genes use []

        bamList = bamFileList[:,0]
        bamsReversed = bamFileList[:,1]
        
        realCov = {}
        
        print ("Calculating real coverage..")
        
        for file_index, bamPath in enumerate(bamList):
            
            baiPath = "".join([bamPath,".bai"])
            
            print (bamPath)
            
            if ( (not os.path.isfile(bamPath)) | (not bamPath[-4:] == ".bam") ):
                print ("Provided path for .bam file is not correct... Please try again...")
                return
        
            if ( not os.path.isfile(baiPath) ):
                
                if ( self.sort ):
                    a = re.search("(\/.*\/)(.*)\.(.*)",bamPath)
                    sortedPath = a[1]+"sorted/"+a[2]+"."+a[3]
                    print ("sorting "+bamPath)
                    if ( not os.path.exists(a[1]+"sorted/") ):
                        os.makedirs(a[1]+"sorted/")
                    ps.sort(bamPath,"-o",sortedPath)
                    bamPath = sortedPath
                
                print ("bai index file not found.. creating index file..")
                ps.index(bamPath)        
                print ("Done.")
            
            bamFile = ps.AlignmentFile(bamPath,"rb")
            genes = self.genes
            
             ## JUST FOR THE STATISTICS ##
            iterations_plus = 0
            iterations_minus = 0
            
            read2 = 0
            read1 = 0
            noread = 0
            duplicates = 0
            unmapped = 0
            secondary = 0
            low_quality_map = 0
            reads = 0
            
            bothStretchCount = 0
            polyTcount = 0
            polyAcount = 0
            polyAdiscarded = 0
            polyTdiscarded = 0
            
            if ( len(self.sampleNames) == 0 ):
                if ( self.extractBarcodes == "fileName" ):
                    sampleName = re.split("/",bamPath)[-1][0:-4]
                    self.sampleNames.append(sampleName) 
                elif ( self.extractBarcodes == "parentDir" ):
                    sampleName = re.split("/",bamPath)[-2]
                    self.sampleNames.append(sampleName)
            else:
                sampleName = self.sampleNames[file_index]
            realCov[sampleName] = {}
            realCov["genePositions"] = {}
            
            for gene in self.genes:
                
                realCov[sampleName][gene] = {}
                realCov[sampleName][gene][self.strands[0]]= {}
                realCov[sampleName][gene][self.strands[1]] = {}

                if (self.UTRmode):
                    gene_annotation = self.selectedAnnotations["whole"][gene]
                else:
                    gene_annotation = self.selectedAnnotations["whole"][gene]

                chromosomeNames = np.unique(gene_annotation["seqname"]).astype(str)
                
                if ( len(chromosomeNames) > 1 ):
                    print ("didn't expect multiple chromosome for one gene..")
                    return
                
                start_genes = np.array(gene_annotation["start"])
                end_genes = np.array(gene_annotation["end"])

                start = np.min(start_genes) - 1
                end = np.max(end_genes)
                                                        
                genePositions = np.arange(start,end)
                
                realCov[sampleName][gene]["+"] = np.zeros(len(genePositions))
                realCov[sampleName][gene]["-"] = np.zeros(len(genePositions))
                
                realCov[sampleName][gene]["numOfReads"] = 0
                
                ##-------------------------##
                for fragmentRead in bamFile.fetch(str(chromosomeNames[0]),start,end):
                    
                    reads += 1
                    
                    ## JUST FOR THE STATISTICS ##
                    
                    if ( fragmentRead.is_secondary ):
                        secondary+=1
                        if ( self.uniqueMap == True ):
                            continue
                    if ( fragmentRead.is_duplicate ):
                        duplicates+=1
                    if ( fragmentRead.is_unmapped ):
                        unmapped+=1
                        if ( self.uniqueMap == True ):
                            continue
                    if ( fragmentRead.mapping_quality < 255 ):
                        low_quality_map += 1
                        if ( self.uniqueMap == True ):
                            continue
                    if ( fragmentRead.is_read1 ):
                        read1+=1
                    elif ( fragmentRead.is_read2 ):
                        read2+=1
                    else:
                        noread+=1
                        
                    ##-------------------------##
                    
                    ## REAL COVERAGE CALCULATION ##
                    
                    relative_indexes = fragmentRead.positions - genePositions[0]
                    
                    last_pos = genePositions[-1] - genePositions[0]
                    ids = np.where(np.logical_and(relative_indexes > 0,relative_indexes <= last_pos))
                    relative_indexes = relative_indexes[ids]
                    
                    if ( self.nonStrandedLib ): ## NON-STRAND-SPECIFIC LIBRARY PREPS

                        ## APPROACH USED FOR OXFORD NANOPORE EXPERIMENT
                        ## principle: 10 consecutive A's within first 50 nt define sense(+) strand, 10 consecutive T's within last 50 nt define antisense(-) strand  

                        readSeq = fragmentRead.query_sequence
                    
                        ab = re.search("T{10,}.*A{10,}",readSeq)
                        if ( not ab == None ):
                            bothStretchCount += 1
                            continue

                        T_match = re.search("T{10,}",readSeq)
                        A_match = re.search("A{10,}",readSeq)
                    
                        if ( ( (T_match == None) & (A_match == None) ) | ( (not T_match == None ) & (not A_match == None) ) ):
                            bothStretchCount += 1
                            continue
                        
                        if ( not T_match == None ):
                            if ( T_match.start() < 50 ): 
                                polyTcount += 1
                                strand = "-"
                            else:
                                polyTdiscarded += 1 
                                continue
                        
                        if ( not A_match == None ):
                            if ( A_match.end() > len(readSeq)-50 ):
                                strand = "+"
                                polyAcount += 1 
                            else:
                                polyAdiscarded += 1 
                                continue

                    else: ## STRAND-SPECIFIC LIBRARY PREPS
                        if (fragmentRead.is_paired): 
                            if (fragmentRead.is_read1):
                                if fragmentRead.is_reverse:
                                    strand = "+"
                                else:
                                    strand = "-"
                            elif (fragmentRead.is_read2):
                                if fragmentRead.is_reverse:
                                    strand = "-"
                                else:
                                    strand = "+"
                        else:
                            if (fragmentRead.is_reverse):
                                strand = "-"
                            else: 
                                strand = "+"
                                
                    realCov[sampleName][gene]["numOfReads"] += 1
                    realCov[sampleName][gene][strand][relative_indexes] += 1
                        
                    length = len(self.sampleNames) - 1

                    ##-------------------------##        
                
                realCov["genePositions"][gene] = genePositions
                
                ## SWAP PLUS AND MINUS STRAND IF FWD LIBRARY PREP
                if ( bamsReversed[file_index] == "REV" ):
                    temp1 = realCov[sampleName][gene]["+"]
                    realCov[sampleName][gene]["+"] = realCov[sampleName][gene]["-"]
                    realCov[sampleName][gene]["-"] = temp1
                    del temp1
                    
                ## COVERAGE FILE CREATION ##
                
                if ( self.createCovFile ):
                    for strand in self.strands:

                        if strand == "-":
                            strandName = "minus" 
                        elif strand == "+":
                            strandName = "plus" 

                        outPath = self.outPath 
                        #outPath += "covs/"

                        if not os.path.exists(outPath):
                            os.makedirs(outPath)

                        file = open("".join([outPath,"gene_coverage_",strandName,".cov"]), "w", encoding="utf-8")

                        for i in np.arange(0,len(realCov[sampleName][gene][strand])):
                            file.write("".join([gene,"\t",str(genePositions[i]),"\t",str(int(realCov[sampleName][gene][strand][i])),"\n"]))

                        file.close()
                
                """
                if (self.createBigWig):
                    print ("TODO: create bigwig file for IGV browser")

                    for strand in self.strands:

                        coverage = realCov[sampleName][gene][strand]

                        start, end, values = self.coverage2intervals(coverage)
                        
                        #with bwig.open(self.outPath + "/igv_coverage/", "w") as bw:
                            #bw.addEntries("chroms", , ends = , values = )
                """

                ##-------------------------##

        if ( self.average ):
            realCov["avg"] = {}
            for gene in self.genes:
                realCov["avg"][gene]={}
                for strand in self.strands:
                    realCov["avg"][gene][strand]={}
                    realCov["avg"][gene][strand] = np.zeros(len(realCov["genePositions"][gene]))
                    for sampleName in self.sampleNames:
                        realCov["avg"][gene][strand] += realCov[sampleName][gene][strand]
                    realCov["avg"][gene][strand] /= len(self.sampleNames)

            self.sampleNames.append("avg")

        self.realCoverage = realCov           
        print ("Done.")

    def calculate_theoretical_coverage(self):   

        print ("Calculating theoretical coverage..") 
        showTransitionLengths = self.showTransitionLengths
        theoretical_coverage = {}  
        
        mix_DF = self.mix

        for mode in list(self.selectedAnnotations):
            
            theoretical_coverage[mode] = {}

            for gene in self.genes:
                
                defaultMolarity = False

                gene_annotation = self.selectedAnnotations[mode][gene] 
        
                maxIndex = np.max(np.array(self.selectedAnnotations["whole"][gene]["end"]))
                minIndex = np.min(np.array(self.selectedAnnotations["whole"][gene]["start"]))
                
                sumCoverages = np.zeros(maxIndex-minIndex+1)
                
                theoretical_coverage[mode][gene] = {}
                theoretical_coverage[mode][gene]["statistics"] = {}
    
                usedTranscriptsAll = np.unique(np.array(mix_DF[1][mix_DF[0] == gene]).flatten())
                
                if ( len(usedTranscriptsAll) == 0 ):
                    if ( self.mixPath!=None ):
                        print ("gene "+gene+" not found in mix file.. considering molarity 1 for all transcripts in gtf file..")
                    
                    defaultMolarity = True
                    usedTranscriptsAll = np.unique(gene_annotation["transcript_id"]).flatten()

                theoretical_coverage[mode][gene]["statistics"]["allUsedTranscriptsCount"] = len(usedTranscriptsAll)
        
                strand_index1 = 0

                maxMolarity = np.zeros(len(self.genes)*2)
                maxCoverage = np.zeros(len(self.genes)*2)

                for direction in self.strands: # looping over strands

                    theoretical_coverage[mode][gene][direction] = {}
                    theoretical_coverage[mode][gene][direction]["transcripts"] = {}

                    # creating array of used genes from mix file

                    usedTranscripts = sorted(gene_annotation["transcript_id"][(gene_annotation["strand"] == direction)&(gene_annotation["transcript_id"].isin(usedTranscriptsAll))].unique())
                    
                    inc = 0
                    if ( len(usedTranscripts) > 0 ):

                        theoretical_coverage[mode][gene][direction]["transcripts"]["names"] = usedTranscripts
                        theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"] = np.zeros((len(usedTranscripts),len(sumCoverages)))
                        
                        for usedTranscript in usedTranscripts:
                            
                            if ( defaultMolarity ):
                                selected_mix_row = np.array([1])
                            else:
                                selected_mix_row = np.array(mix_DF[2][mix_DF[1]  ==  usedTranscript]).flatten()
                                
                            if ( len(selected_mix_row) > 1 ):
                                print ("Multiple rows for transcript "+usedTranscript+" in mix file!! Taken first row value..")
                                
                            starts = np.array(gene_annotation["start"][ (gene_annotation["transcript_id"]  ==  usedTranscript) & (gene_annotation["strand"]  ==  direction)] ).flatten()
                            starts -= minIndex 
                            ends = np.array(gene_annotation["end"][( gene_annotation["transcript_id"]  ==  usedTranscript) & (gene_annotation["strand"]  ==  direction) ]).flatten()
                            ends -= minIndex + 1
                            molarity = float(selected_mix_row)
                            transcriptIndex = int(usedTranscripts.index(usedTranscript))

                            for i in np.arange(0,len(starts)):
                                theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"][transcriptIndex,starts[i]:ends[i]] += molarity # increasing coverage weighted by the molarity

                                if ( (showTransitionLengths) & (i == 0) ):
                                    exonTotalLength = np.sum(ends - starts)
                                    startDistance,_,_ = self.calc_start_transition(exonTotalLength,30,25)
                                    startDistance = np.flip(startDistance, axis = 0)[0:(ends[i]-starts[i])]
                                    theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"][transcriptIndex,starts[i]:ends[i]] *= startDistance

                                if ( (showTransitionLengths) & (i == len(starts)-1) ):
                                    exonTotalLength = np.sum(ends - starts)
                                    startDistance,_,_ = self.calc_start_transition(exonTotalLength,25,25)
                                    startDistance = startDistance[len(startDistance)-(ends[i]-starts[i]):len(startDistance)+1]
                                    theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"][transcriptIndex,starts[i]:ends[i]] *= startDistance

                            inc += 1

                        theoretical_coverage[mode][gene][direction]["transcripts"]["starts"] = starts
                        theoretical_coverage[mode][gene][direction]["transcripts"]["ends"] = ends
                        theoretical_coverage[mode][gene][direction]["coverage"] = np.sum(theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"],axis=0) 
                                
                        sumCoverages += theoretical_coverage[mode][gene][direction]["coverage"]
                        
                        maxMolarity[strand_index1] = np.max(theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"])
                        maxCoverage[strand_index1] = np.max(theoretical_coverage[mode][gene][direction]["coverage"])
                    else: 
                        theoretical_coverage[mode][gene][direction]["transcripts"]["names"] = []
                        theoretical_coverage[mode][gene][direction]["transcripts"]["matrix"] = []
                        theoretical_coverage[mode][gene][direction]["coverage"] = []
                        theoretical_coverage[mode][gene][direction]["UTR_coverage"] = []

                    strand_index1 += 1

                theoretical_coverage[mode][gene]["statistics"]["maxMolarity"] = np.max(maxMolarity)
                theoretical_coverage[mode][gene]["statistics"]["maxCoverage"] = np.max(maxCoverage)

                starts,ends = self.returnChanges_greater0(sumCoverages)
                
                theoretical_coverage[mode][gene]["segments"]={}
                theoretical_coverage[mode][gene]["segments"]["starts"] = starts
                theoretical_coverage[mode][gene]["segments"]["ends"] = ends
        
        #if (self.UTRmode):
        #    for gene in self.genes:
        #        theoretical_coverage["UTR_only"][gene]["segments"]["ends"] = theoretical_coverage["whole"][gene]["segments"]["ends"]
        #        theoretical_coverage["UTR_only"][gene]["segments"]["starts"] = theoretical_coverage["whole"][gene]["segments"]["starts"]

        self.theoretical_coverage = theoretical_coverage

        print ("Done.")

    def calculate_stats(self):
        
        stats = {}
        stats["COD"] = {}
        stats["scaleCoef"] = {}
        
        if ( self.CODCalc ):
            for gene in self.genes:
                
                stats["COD"][gene] = {}
                stats["scaleCoef"][gene] = {}

                for sampleName in self.sampleNames:

                    stats["COD"][gene][sampleName] = {}
                    stats["scaleCoef"][gene][sampleName] = {}
                    
                    for strand in self.strands:

                        if (self.UTRmode):
                            mode = 'UTR_only'
                        else:
                            mode = 'whole'

                        coverage = self.theoretical_coverage[mode][gene][strand]["coverage"]
                        real_coverage = self.realCoverage[sampleName][gene][strand]

                        if ( len(coverage) != 0 ):

                            ## if UTR mode, CoD only within 3' UTR region, the rest is ignored 
                            if ( mode == "UTR_only" ):
                                real_coverage = np.array(coverage > 0).astype(float) * real_coverage

                            COD, scaling_coefficient = self.CoD(real_coverage,coverage)
                            stats["COD"][gene][sampleName][strand] = COD
                            stats["scaleCoef"][gene][sampleName][strand] = scaling_coefficient 
                        else:
                            stats["COD"][gene][sampleName][strand] = "NaN"
                            stats["scaleCoef"][gene][sampleName][strand] = "NaN"

            ## WRITING COD TABLE ##
            
            outPath = self.outPath

            if ( not os.path.exists(outPath) ):
                os.makedirs(outPath)

            with open(outPath+"gene_"+gene+"_CoD_table.csv", 'w') as COD_file:
                header = ";"
                for sampleName in self.sampleNames:
                    for strand in self.strands:
                        header += sampleName + " (" + strand + ");"

                COD_file.write(header+"\n")
 
                for gene in self.genes:
                    row = gene + ";"
                    for sampleName in self.sampleNames:
                        for strand in self.strands:
                            row += str(stats["COD"][gene][sampleName][strand]) + ";"
                        row += "\n"
                    COD_file.write(row)

            self.stats = stats
    
    def CoD(self,realCov,theoryCov):
        ## This function calculate CoD from 2 vectors (real and expected coverage)
        ## and returns coefficient for scaling and CoD
        ## input: 2x ndarray
        ## output: CoD and scaling coefficient 
        
        sumBAC_theory = np.sum(theoryCov) # sumBAC (sum of base annotation count), sum of expected coverage
        
        sumCoverage_real = np.sum(realCov) # sum of real coverage

        scaling_coefficient = sumCoverage_real/sumBAC_theory # calculation of scaling factor for real coverage (for comparison)
        
        ## Correction of 0 factor
        if ( scaling_coefficient == 0 ): 
            scaling_coefficient = 1
            
        ## Scaling
        coverageReal_scaled = realCov/scaling_coefficient
            
        ## COD computation 
        difference = theoryCov - coverageReal_scaled
        COD = np.sum((difference)**2)/sumBAC_theory
        COD = np.round(COD,4)
        
        return COD, scaling_coefficient

    def coverage2intervals(self,coverage):
        
        coverage = np.append(0,coverage)
        coverage = np.append(coverage,0)

        differences = np.diff(coverage)
        start = np.argwhere(differences > 0).flatten()
        end = np.argwhere(differences < 0).flatten()

        values = coverage[start+1]

        return start, end, values
        
    def returnChanges_greater0(self,vector):
        ## This function gives indexes of input vector where value differs from 0
        ## input: ndarray num vector
        ## output: nddarray of position of positive or negative differences
        
        # extending vector with zeros to determine even start and ends # 
        vector=np.append(0,vector)
        vector=np.append(vector,0) 
        # getting position of positive and negative differences # 
        differences = np.diff((vector>0).astype(int))
        start = np.argwhere(differences == 1).flatten()
        end = np.argwhere(differences == -1).flatten()
        
        return start,end

    def returnChanges_all(self,vector):
        ## This function gives indexes of non-zero differences
        ## input: ndarray num vector
        ## output: nddarray of position of differences
        differences = np.argwhere((np.diff(vector)!=0).astype(int)).flatten()    
        differences = np.append(0,differences)
        differences = np.append(differences,len(vector))
        
        return differences
    
    def calc_start_transition(self,transcriptLength,mean=200,std=80):
        ## This function estimates probability distribution of a fragment according to its length,
        ## input: length offragment, mean and standard deviation is inserted
        ## output: transition length according to probability, gaussian distribution of length probability
        ##         and length where the prob. is high enough
        
        lengthProbability = norm.pdf(np.arange(1,transcriptLength+1),mean,std)
        
        lengthProbability /= np.sum(lengthProbability)
        
        fragmentProbability = np.zeros((transcriptLength,transcriptLength))
        
        for length in np.arange(1,transcriptLength):
            fragmentProbability[length,0:transcriptLength-length+1] = lengthProbability[length]/(transcriptLength-length+1)
        
        startDistance = np.sum(fragmentProbability,axis=0)
        
        adjustedLength = np.sum(lengthProbability*(transcriptLength-np.arange(1,transcriptLength+1)+1))
        
        startDistance /= np.max(startDistance)
        
        return startDistance, lengthProbability, adjustedLength


    def draw_header(self):
        print ("")

    def draw_transcript_chart(self):
        print ("")

    def draw_coverage_chart(self):
        print ("")

    def draw_text(self, text = None, x = None, y = None, font_size = 28, rotation = 0, color_rgb = (0,0,0), h_align = "center", v_align = "center"):
        
        ctx = self.cr
        ctx.set_source_rgb(color_rgb[0],color_rgb[1],color_rgb[2])
        ctx.set_font_size(font_size)
        
        (_, _, width, height, _, _) = ctx.text_extents(text) 

        if h_align == "center":
            x = x + width/2 
        elif h_align == "right":
            x = x + width

        if v_align == "center":
            y = y + height/2
        elif v_align == "bottom":
            y = y + height
        elif v_align == "top":
            y = y

        ctx.move_to(x,y)
        cr.show_text(text)


        
    def visualize2png(self,relative_scaling,colorSpec=[]):
    
        print ("Creating visualizations...")
        
        showReferenceSegments = False
        coverage_height_perc = 40
        coverage_width_perc = 85
            
        width_paper = 1980
        height_paper = 1224
        
        formats = ["png"]
        
        TITLE = ""
        EXPERIMENT_NAME = self.experimentName
        NAME = "Sample:  "
        TITLE = "Gene:  " 
        SUBTITLE = "Experiment:  "  

        if ( len(colorSpec) == 0 ):
            colorsReal = [[0.4, 0.4, 0.4, 0.7],
                        [0, 1, 0, 0.6],
                        [1,0,0,0.6],
                        [1,1,0.4,1],
                        [0,1,1,1],
                        [0.8,0.4,0.4,1],
                        [1,0.5,0.4,1],
                        [1,0,1,1]]
        else:
            colorsReal = np.array([colors.to_rgba(colorSpec)])
        
        for gene in self.genes:
            
            #surfaces = []
            
            #outPath = self.outPath+"coverage_plots/"+self.sampleNames[0]+"/"+gene+"/"
            #outPath = self.outPath

            if ( not os.path.exists(outPath) ):
                os.makedirs(outPath)

            for formatIndex in np.arange(0,len(formats)):
                
                """

                if ( formats[formatIndex] == "svg" ):
                    surfaces.append(cairo.SVGSurface("".join([outPath,"coverage_",NAME,"_",gene,".svg"]), width_paper, height_paper))
                elif ( formats[formatIndex] == "pdf" ):
                    surfaces.append(cairo.PDFSurface("".join([outPath,"coverage_",NAME,"_",gene,".pdf"]), width_paper, height_paper))
                elif ( formats[formatIndex] == "png" ):
                    surfaces.append(cairo.ImageSurface(cairo.FORMAT_ARGB32, width_paper, height_paper))
                elif ( formats[formatIndex] == "ps" ):
                    surfaces.append(cairo.PSSurface("".join([outPath,"coverage_",NAME,"_",gene,".ps"]), width_paper, height_paper))
    
                cr = cairo.Context(surfaces[formatIndex])

                """
                
                cr.set_source_rgb(1, 1, 1)
                
                modes = np.flip(np.sort(list(self.selectedAnnotations)), axis = 0)

                cr.rectangle(0,0,width_paper,height_paper)
                cr.fill()
                perc_height=15
                page_half_x = width_paper/4

                ## PREPARE FUNCTION FOR CALCULATING COVERAGE SCALING

                if ( EXPERIMENT_NAME != "" ):
                        cr.set_font_size(28)
                        (_, _, width, height, _, _) = cr.text_extents(SUBTITLE) 
                        cr.set_source_rgb(0.2,0.2,0.2)
                        cr.move_to(page_half_x-width,100)
                        SUBTITLE += EXPERIMENT_NAME
                        cr.show_text(SUBTITLE)    

                cr.set_font_size(28)
                cr.set_source_rgb(0.2,0.2,0.2)
                namePositionY = height_paper/12*perc_height/100
                cr.move_to(page_half_x - width,60)
                cr.show_text(NAME)
                NAME2 = "  " + self.sampleNames[0]
                cr.move_to(page_half_x + 4,60)
                cr.show_text(NAME2)


                cr.select_font_face("Sans", cairo.FONT_SLANT_NORMAL,
                cairo.FONT_WEIGHT_NORMAL)
                cr.set_font_size(28)
                cr.set_source_rgb(0.2,0.2,0.2)
                cr.move_to(page_half_x - width,140)
                cr.show_text(TITLE)

                TITLE2 = "  " + gene 
                if self.UTRmode:
                    TITLE2 += " - 3' UTR coverage" 
                else:
                    TITLE2 += " - whole coverage"

                cr.move_to(page_half_x + 4, 140)
                cr.show_text(TITLE2)

                #cov_value_to_draw = (lambda coverage_value,scalingCoef: )

                for mode in modes:

                    if ( self.UTRmode and mode == "whole"):
                        gray_background = True
                    else:
                        gray_background = False
                
                    theoretical_coverage = self.theoretical_coverage
                    coverage_max = theoretical_coverage[mode][gene]["statistics"]["maxCoverage"] 
                    
                    """
                    sampleText = "Sample: " + self.sampleNames[0]
                    (_, _, width, height, _, _) = cr.text_extents(sampleText)
                    cr.set_source_rgb(0.2,0.2,0.2)
                    cr.move_to(width_paper/2-width/2,120)
                    cr.show_text(sampleText)
                    """                    

                    coveragePlot_width = coverage_width_perc*width_paper/100
                    coveragePlot_offsetX = (width_paper-coveragePlot_width)/2
                    coveragePlot_height = coverage_height_perc*height_paper/100
                    coveragePlot_offsetY = 150

                    coveragePlot2_Y = coveragePlot_height+coveragePlot_offsetY+220 #height_paper/5
                    coveragePlot2_height = 170

                    offsetX_text = 20
                    offsetY_lines = coveragePlot_height/10
                    yLine = coveragePlot_offsetY+offsetY_lines

                    axis_description_x = coveragePlot_offsetX + coveragePlot_width - 20
                    axis_description_y = coveragePlot2_Y - coveragePlot2_height/6 

                    exonHeight = offsetY_lines*0.4

                    cr.set_source_rgb(0, 0, 0)
                    cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL,
                            cairo.FONT_WEIGHT_NORMAL)

                    allTranscriptsCount = theoretical_coverage[mode][gene]["statistics"]["allUsedTranscriptsCount"]
                    maxMolarity = theoretical_coverage[mode][gene]["statistics"]["maxMolarity"]
                    
                    StatisticsY = coveragePlot2_Y - coveragePlot2_height / 2
                    StatisticsX = coveragePlot_offsetX / 4


                    segmentStarts = theoretical_coverage["whole"][gene]["segments"]["starts"]
                    segmentEnds = theoretical_coverage["whole"][gene]["segments"]["ends"]
                    segmentCount = len(segmentStarts)
                    exonLengths = segmentEnds-segmentStarts

                    for strand in ["-","+"]:

                        matrix2draw = theoretical_coverage[mode][gene][strand]["transcripts"]["matrix"]

                        if ( len(matrix2draw) == 0 ):
                            continue

                        sirvTranscripts = theoretical_coverage[mode][gene][strand]["transcripts"]["names"]
                        starts = theoretical_coverage["whole"][gene][strand]["transcripts"]["starts"]
                        ends = theoretical_coverage["whole"][gene][strand]["transcripts"]["ends"]
                        
                        cr.set_font_size(exonHeight*12/10)
                        (x, y, width, height, dx, dy) = cr.text_extents(sirvTranscripts[0]) 
                        startX_text = coveragePlot_offsetX+offsetX_text
                        startX_line = startX_text+offsetX_text+width
                        endX_line = coveragePlot_offsetX+coveragePlot_width-offsetX_text

                        if (gray_background):
                            cr.set_source_rgba(0,0,0,0.5)
                        else:
                            cr.set_source_rgb(0,0,0)

                        cr.set_line_cap(cairo.LINE_CAP_ROUND)

                        lineInnerHeight = coveragePlot_height-offsetY_lines
                        lineLength = endX_line-startX_line

                        for i in np.arange(0,len(sirvTranscripts)): 

                            
                            cr.set_font_size(exonHeight*12/10)
                            (_, _, width, height, _, _) = cr.text_extents(sirvTranscripts[i]) 

                            cr.move_to(startX_text,yLine+height/2)
                            
                            if (gray_background):
                                cr.set_source_rgba(0,0,0,0.5)
                            else:
                                cr.set_source_rgb(0,0,0)
                            
                            cr.show_text(sirvTranscripts[i])
                            
                            if (not gray_background):
                                cr.set_line_width(2)
                                cr.move_to(startX_line,yLine)
                                cr.line_to(endX_line,yLine)
                                cr.stroke()

                            segmentGaps = lineLength/50
                            startSegmentsDraw = startX_line + segmentGaps
                            sumOfLengths = np.sum(exonLengths)
                            segmentExonsLength_draw = (lineLength)-(segmentCount+1)*segmentGaps

                            for j in np.arange(0,segmentCount):

                                segmentIndexStart = theoretical_coverage["whole"][gene]["segments"]["starts"][j]
                                segmentIndexEnd = theoretical_coverage["whole"][gene]["segments"]["ends"][j]

                                exonSegment = matrix2draw[i,segmentIndexStart:segmentIndexEnd]

                                exonHeight_proportional = np.max(matrix2draw[i,:])/maxMolarity*exonHeight   

                                basic_scaling = (coveragePlot2_height)/(coverage_max*exonHeight)

                                exonStarts, exonEnds = self.returnChanges_greater0(exonSegment)

                                segmentLength = exonLengths[j]/sumOfLengths*(segmentExonsLength_draw)

                                if ( showReferenceSegments ):
                                    
                                    if (gray_background):
                                        cr.set_source_rgba(0.95,0.95,0.95,0.5)
                                    else:
                                        cr.set_source_rgba(0.95,0.95,0.95)
                                    
                                    cr.rectangle(startSegmentsDraw,yLine-exonHeight_proportional/2,segmentLength,exonHeight_proportional)
                                    cr.fill_preserve()

                                    if (gray_background):
                                        cr.set_source_rgba(0.2,0.2,0.2,0.5)
                                    else:
                                        cr.set_source_rgb(0.2,0.2,0.2)
                            
                                    cr.set_line_width(1)
                                    cr.stroke()

                                if ( len(exonStarts) > 0 ):
                                    if ( strand == "+" ):
                                        if (gray_background):
                                            cr.set_source_rgba(0.3,0.8,0.3,0.5)
                                        else:
                                            cr.set_source_rgb(0.3,0.8,0.3)
                                        
                                    elif ( strand == "-" ):
                                        if (gray_background):
                                            cr.set_source_rgba(0.3,0.7,1,0.5)
                                        else:
                                            cr.set_source_rgb(0.3,0.7,1)

                                    for exonIndex in range(0,len(exonStarts)):
                                        exonStartPosX = startSegmentsDraw+exonStarts[exonIndex]/len(exonSegment)*segmentLength
                                        exonLength = startSegmentsDraw+exonEnds[exonIndex]/len(exonSegment)*segmentLength-exonStartPosX
                                        cr.rectangle(exonStartPosX,yLine-exonHeight_proportional/2,exonLength,exonHeight_proportional)# Rectangle(x0, y0, x1, y1)                    cr.stroke()

                                    cr.fill_preserve()

                                    if (gray_background):
                                        cr.set_source_rgba(0.2,0.2,0.2,0.5)
                                    else:
                                        cr.set_source_rgb(0.2,0.2,0.2)

                                    cr.set_line_width(0.04)
                                    cr.stroke()
                                
                                ## DURING FIRST SEGMENT OF TRANSCRIPT ANNOTATION, PROJECT COVERAGES TO THE LOWER COVERAGE CHART

                                if ( ( self.UTRmode and mode == "UTR_only" ) or ( self.UTRmode == False and mode == "whole" ) ):

                                    if ( i == 0 ):
                                        
                                        # draw x-axis legend
                                        noNullStrand = [k for k in ["+","-"] if len(theoretical_coverage[mode][gene][k]["transcripts"]["names"])>0][0]
                                        if strand == noNullStrand:

                                            cr.set_font_size(16)
                                            cr.set_source_rgba(0, 0, 0, 1)
                                            cr.set_line_width(4)

                                            x_axis_coords = str(self.realCoverage["genePositions"][gene][segmentStarts[j]])
                                            (_, _, _, height_text, _, _) = cr.text_extents(x_axis_coords) 
                                            cr.move_to(startSegmentsDraw - height_text/2,coveragePlot2_Y+coveragePlot2_height+110)
                                            cr.rotate(90*np.pi/180)
                                            cr.show_text(x_axis_coords)
                                            cr.rotate(-90*np.pi/180)
                                            
                                            x_axis_coords = str(self.realCoverage["genePositions"][gene][segmentEnds[j]])
                                            (_, _, _, height_text, _, _) = cr.text_extents(x_axis_coords) 
                                            cr.move_to(startSegmentsDraw+segmentLength - height_text/2,coveragePlot2_Y+coveragePlot2_height+110)
                                            cr.rotate(90*np.pi/180)
                                            cr.show_text(x_axis_coords)
                                            cr.rotate(-90*np.pi/180)

                                            cr.set_source_rgba(0.2, 0.2, 0.2, 0.7)

                                            cr.move_to(startSegmentsDraw,coveragePlot2_Y+coveragePlot2_height+100)
                                            cr.line_to(startSegmentsDraw+segmentLength,coveragePlot2_Y+coveragePlot2_height+100)
                                            cr.move_to(startSegmentsDraw,coveragePlot2_Y+coveragePlot2_height+105)
                                            cr.line_to(startSegmentsDraw,coveragePlot2_Y+coveragePlot2_height+95)
                                            cr.move_to(startSegmentsDraw+segmentLength,coveragePlot2_Y+coveragePlot2_height+105)
                                            cr.line_to(startSegmentsDraw+segmentLength,coveragePlot2_Y+coveragePlot2_height+95)
                                            
                                            cr.stroke()

                                        coverageSegment = theoretical_coverage[mode][gene][strand]["coverage"][segmentIndexStart:segmentIndexEnd]
                                        coverageSegment = np.append(0,coverageSegment)
                                        coverageSegment = np.append(coverageSegment,0)
                                        indexes4change = self.returnChanges_all(coverageSegment) 

                                        cr.move_to(startSegmentsDraw,coveragePlot2_Y)

                                        posits = 0
                                        legend_positionX = 20
                                        legend_positionY = 600
                                        
                                        if ( self.average ):
                                            sampleNames=["avg"]
                                        else:
                                            sampleNames=self.sampleNames
                                        
                                        for sampleName in sampleNames:

                                            segment_start = theoretical_coverage["whole"][gene]["segments"]["starts"][j]
                                            segment_end = theoretical_coverage["whole"][gene]["segments"]["ends"][j]

                                            
                                            realCoverageSegment = self.realCoverage[sampleName][gene][strand][segment_start:segment_end]

                                            cr.move_to(startSegmentsDraw,coveragePlot2_Y)

                                            limit_draw_lower = coveragePlot2_Y + coveragePlot2_height + 50
                                            limit_draw_upper = coveragePlot2_Y - coveragePlot2_height - 50
                                            
                                            scaling_coefficients = self.stats["scaleCoef"][gene][sampleName]

                                            if ('NaN' not in scaling_coefficients.values()):
                                                coverage_scaling_coef = max(list(scaling_coefficients.values()))
                                            elif (scaling_coefficients['+'] == 'NaN'):
                                                coverage_scaling_coef = scaling_coefficients['-']
                                            elif (scaling_coefficients['-'] == 'NaN'):
                                                coverage_scaling_coef = scaling_coefficients['+']

                                            max_visible_read = {}
                                            if (scaling_coefficients["+"] == "NaN"):
                                                max_visible_read["+"] = None
                                            else:
                                                max_visible_read["+"] = abs(coveragePlot2_Y - limit_draw_upper) * coverage_scaling_coef * relative_scaling / basic_scaling / exonHeight
                                            
                                            if (scaling_coefficients["-"] == "NaN"):
                                                max_visible_read["-"] = None
                                            else:
                                                max_visible_read["-"] = abs(coveragePlot2_Y - limit_draw_lower) * coverage_scaling_coef * relative_scaling / basic_scaling / exonHeight
                                            
                                            realCoverageSegment = np.append(0,realCoverageSegment)
                                            
                                            for realCovarageIndex in np.arange(0,len(realCoverageSegment)):
                                                if ( strand == "-" ):
                                                    value = -realCoverageSegment[realCovarageIndex]
                                                    color = [0.2,0.2,0.2,0.5]

                                                elif ( strand == "+" ):
                                                    value = realCoverageSegment[realCovarageIndex]
                                                    color = [0.2,0.2,0.2,0.5]

                                                draw_X = startSegmentsDraw + realCovarageIndex/len(realCoverageSegment)*segmentLength
                                                draw_Y = coveragePlot2_Y - value / coverage_scaling_coef * relative_scaling * basic_scaling * exonHeight
                                                
                                                if (draw_Y > limit_draw_lower):
                                                    draw_Y = limit_draw_lower

                                                if (draw_Y < limit_draw_upper):
                                                    draw_Y = limit_draw_upper

                                                cr.line_to(draw_X,draw_Y)

                                            cr.line_to(startSegmentsDraw+segmentLength,coveragePlot2_Y)

                                            cr.close_path()
                                            color = colorsReal[0]
                                            cr.set_source_rgba(color[0],color[1],color[2],color[3])
                                            cr.fill_preserve()

                                            cr.set_source_rgb(colorsReal[posits][0],colorsReal[posits][1],colorsReal[posits][2])
                                            cr.set_line_width(2)
                                            cr.stroke()

                                            

                                            ## DRAW EXCEEDED CUTS

                                            scaled_coverage = realCoverageSegment / ( coverage_scaling_coef * relative_scaling ) * basic_scaling * exonHeight
                                            max_possible_coverage = limit_draw_upper / ( coverage_scaling_coef * relative_scaling ) * basic_scaling * exonHeight
                                            min_possible_coverage = limit_draw_lower / ( coverage_scaling_coef * relative_scaling ) * basic_scaling * exonHeight
                                            
                                            triangle_shift = 5
                                            triangle_height = 10

                                            if (strand == '-'):
                                                exceedLimitSegment = coveragePlot2_Y + scaled_coverage
                                                exceedLimitSegment = exceedLimitSegment > limit_draw_lower

                                                limit_y = limit_draw_lower

                                            elif (strand == '+'):
                                                exceedLimitSegment = coveragePlot2_Y - scaled_coverage
                                                exceedLimitSegment = exceedLimitSegment < limit_draw_upper

                                                limit_y = limit_draw_upper
                                                triangle_shift = -triangle_shift
                                                triangle_height = -triangle_height

                                            start_exceed, end_exceed = self.returnChanges_greater0(exceedLimitSegment)
                                            
                                            if (len(start_exceed) > 0):

                                                # merge exceeding segments if they are close together to prevent unclear overlap of marks
                                                
                                                distances = start_exceed[1:] - end_exceed[0:-1]

                                                start_exceed_mark = list(start_exceed)
                                                end_exceed_mark = list(end_exceed)
                                                removed = 0 

                                                for dist_index, dist in enumerate(distances):
                                                    if (dist < 20):
                                                        del start_exceed_mark[dist_index + 1 - removed]
                                                        del end_exceed_mark[dist_index - removed]
                                                        removed += 1
                                                
                                                cr.set_source_rgba(1, 0.2, 0.2, 0.7)
                                                cr.set_line_width(2)                                                

                                                # draw red line for exceeding coverage
                                                for index in range(len(start_exceed)):
                                                    start = startSegmentsDraw + start_exceed[index] / len(realCoverageSegment) * segmentLength
                                                    end = startSegmentsDraw + end_exceed[index] / len(realCoverageSegment) * segmentLength

                                                    cr.move_to(start, limit_y)
                                                    cr.line_to(end, limit_y)
                                                    cr.stroke()

                                                # draw mark for exceeding coverage

                                                for index in range(len(start_exceed_mark)):
                                                    start = startSegmentsDraw + start_exceed_mark[index] / len(realCoverageSegment) * segmentLength
                                                    end = startSegmentsDraw + end_exceed_mark[index] / len(realCoverageSegment) * segmentLength
                                                    middle = ( end + start ) / 2
                                                    
                                                    half_width = 5    
                                                    cr.move_to(middle - half_width, limit_y + triangle_shift)
                                                    cr.line_to(middle + half_width, limit_y + triangle_shift)
                                                    cr.line_to(middle, limit_y + triangle_shift + triangle_height)
                                                    cr.close_path()
                                                    cr.stroke()
                                                
                                            cr.set_source_rgba(colorsReal[posits][0],colorsReal[posits][1],colorsReal[posits][2])
                                            cr.move_to(legend_positionX,legend_positionY)
                                            cr.set_font_size(20)
                                            
                                            #if ( self.average ):
                                            #    text = ""
                                            #    for name in self.sampleNames:
                                            #        text += name
                                            #        text += " "
                                                    
                                            #    cr.show_text(text)
                                            #else:
                                            #    cr.show_text(sampleName)
                                            
                                            posits+=1
                                            legend_positionY += 40

                                            cr.move_to(startSegmentsDraw,coveragePlot2_Y)

                                        for k in np.arange(0,len(indexes4change)-1):

                                            if ( strand == "-" ):
                                                value = -coverageSegment[indexes4change[k]+1] 
                                                color = [0,0.5,1,0.2]
                                            elif ( strand == "+" ):
                                                value = coverageSegment[indexes4change[k]+1]
                                                color = [0,0.7,0,0.2]

                                            draw_startX = startSegmentsDraw+indexes4change[k]/len(coverageSegment)*segmentLength
                                            draw_startX_next = startSegmentsDraw+indexes4change[k+1]/len(coverageSegment)*segmentLength
                                            draw_startY = coveragePlot2_Y-value/coverage_max*coveragePlot2_height

                                            cr.line_to(draw_startX,draw_startY)
                                            cr.line_to(draw_startX_next,draw_startY)

                                        cr.line_to(startSegmentsDraw+segmentLength,coveragePlot2_Y)
                                        cr.close_path()

                                        cr.set_source_rgba(color[0],color[1],color[2],color[3])
                                        cr.fill_preserve()

                                        cr.set_source_rgb(0,0,0)
                                        cr.set_line_width(1)
                                        cr.stroke()

                                startSegmentsDraw += segmentLength + segmentGaps   

                            yLine += (lineInnerHeight)/allTranscriptsCount

                # DRAW Y AXIS TEXT AND ARROWS 
                
                offset_y_arrow = 6

                if (max_visible_read["+"] is not None):
                    cr.move_to(startX_line - 10, limit_draw_upper )
                    cr.line_to(startX_line - 10, coveragePlot2_Y - offset_y_arrow )
                    cr.move_to(startX_line - 10 + 4, limit_draw_upper + offset_y_arrow)
                    cr.line_to(startX_line - 10, limit_draw_upper)
                    cr.line_to(startX_line - 10 - 4, limit_draw_upper + offset_y_arrow)

                    y_axis_label_plus = str( int(np.round(max_visible_read["+"])) ) + " reads"
                    (_, _, width_text, height_text, _, _) = cr.text_extents(y_axis_label_plus) 
                    y_label_position_plus = coveragePlot2_Y - (coveragePlot2_Y - limit_draw_upper)/2 - width_text/2
                    cr.move_to(startX_line - 35, y_label_position_plus)
                    cr.rotate(90*np.pi/180)
                    cr.show_text(y_axis_label_plus)
                    cr.rotate(-90*np.pi/180)

                if (max_visible_read["-"] is not None):
                    cr.move_to(startX_line - 10, coveragePlot2_Y + offset_y_arrow )
                    cr.line_to(startX_line - 10, limit_draw_lower )
                    cr.move_to(startX_line - 10 + 4, limit_draw_lower - offset_y_arrow)
                    cr.line_to(startX_line - 10, limit_draw_lower)
                    cr.line_to(startX_line - 10 - 4, limit_draw_lower - offset_y_arrow)

                    y_axis_label_minus = str( int(np.round(max_visible_read["-"])) ) + " reads"
                    (_, _, width_text, height_text, _, _) = cr.text_extents(y_axis_label_minus) 
                    y_label_position_minus = (limit_draw_lower - coveragePlot2_Y)/2 + coveragePlot2_Y - width_text/2
                    cr.move_to(startX_line - 35, y_label_position_minus)
                    cr.rotate(90*np.pi/180)
                    cr.show_text(y_axis_label_minus)
                    cr.rotate(-90*np.pi/180)

                cr.set_source_rgba(0.2, 0.2, 0.2, 0.7)
                cr.set_line_width(4)
                cr.stroke()

                cr.set_source_rgba(0.2, 0.2, 0.2, 0.7)
                cr.set_font_size(24)       
                
                # DRAW COD VALUES
                for sampleName in self.sampleNames:
                    for strand in ["+","-"]:
                        COD = self.stats["COD"][gene][self.sampleNames[0]][strand]
                        if ( COD is not "NaN" ):
                            COD = str(COD)

                            if ( strand == "-" ):
                                color = [0,0.5,1]
                            elif ( strand == "+" ):
                                color = [0,0.7,0]

                            cr.set_source_rgba(color[0],color[1],color[2])
                            cr.move_to(StatisticsX,StatisticsY)
                            cr.set_font_size(24)
                            cr.show_text("".join(["CoD (",strand,"): ",COD]))

                        StatisticsY += coveragePlot_height / 2

                # DRAW COVERAGE AXIS STRANDNESS INDICATOR
                for strand in ["+","-"]:
                    if (max_visible_read[strand] is not None):
                        cr.set_source_rgba(0.2,0.2,0.2,1)
                        cr.move_to(axis_description_x, axis_description_y)
                        cr.set_font_size(22)
                        if ( strand == "+" ):
                            axis_description_text = "sense"
                        elif ( strand == "-" ):
                            axis_description_text = "antisense"

                        cr.show_text(axis_description_text + " ("+ strand +")")
                    axis_description_y += coveragePlot_height/6 + 22
            
                cr.show_page()

                if ( formats[formatIndex]=="png" ):
                    surfaces[formatIndex].write_to_png("".join([outPath,"gene_",gene,"_coverage_plot.png"]))
                
                surfaces[formatIndex].finish()

            del surfaces

        print ("Coverage files can be found in: " + self.outPath)
        print ("Done.")
        return