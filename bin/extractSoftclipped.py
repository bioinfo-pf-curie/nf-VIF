#!/usr/bin/env python

# Author(s): Nicolas Servant and Tina Alaeitabar
# Contact: nicolas.servant@curie.fr

# Copyright Institut Curie 2019
# This software is a computer program whose purpose is to analyze high-throughput sequencing data.
# You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
# The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
# Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
# The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

""" 
Scripts to manage soft-clipped reads in BAM file
"""
import argparse
import sys
import os
import re
import pysam
import pandas as pd

def isHardClipped(read):
    """Check if a read is hard-clipped"""
    pattern = re.compiler('(\w+)H')
    if pattern.match(read.cigarstring):
        return True
    else:
        return False

def isSoftClipped(read):
    """Check if a read is soft-clipped"""
    pattern = re.compile('(\w+)S')
    if pattern.match(read.cigarstring):
        return True
        
    else:
        return False

def ReverseComplement(Pattern):
    """Reverse Complement of a sequence"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(Pattern))
    return reverse_complement


def getSoftClippedCoord(read, annot=False):
    """Get the clipped position on the read from the CIGAR string
       annot : if true, distinguish 5'/3' end coordinates
       5prime = read clipped at the beginning (5' end)
       3prime = read clipped in the end (3' end)
    """
    idx = 0
    pos = []
    label = []
    ## 0 - Match / 1 - Insertion / 2 - Deletion / 3 - Skipped / 4- Soft clipped / 5- Hard clipped / 6- Padding
    for (cigarType,cigarLength) in read.cigartuples:
        try:
            if(cigarType == 4):
                pos.append((idx, idx+cigarLength))
                if annot:
                    if idx == 0:
                        label.append('5prime')
                    elif idx+cigarLength == len(read.seq):
                        label.append('3prime')
                    else:
                        label.append(None)
                else:
                    label.append(None)
                        
            if(cigarType != 2 and cigarType != 5):
                idx += cigarLength
        except:
            sys.stderr.write("Error in CIGAR parsing for read", read.name)
    return (pos, label)
    

def getGenomicBkp(read, excludeBorders=False, reflen=None):
    """Return the breakpoint(s) position(s) on the reference (HPV) genome 
       excludeBorders : remove extreme values at 0 and len(ref) position
    """
    idx = read.reference_start
    pos = []

    if reflen is None and excludeBorders:
        sys.stderr.write("Warning: Reference length is not defined. Cannot exclude borders !")
        excludeBorders = False

    for (cigarType,cigarLength) in read.cigartuples:
        try:
            if(cigarType == 4):
                if idx == 0 and excludeBorders:
                    continue
                elif idx == reflen and excludeBorders:
                    continue
                else:
                    pos.append(idx)
                        
            elif(cigarType != 2 and cigarType != 5):
                idx += cigarLength
        except:
            sys.stderr.write("Error in CIGAR parsing for read", read.name)
    if len(pos) > 0:
        return pos
    else:
        return None
    

def getClippedSeq(read, coords, revert=False):
    """Return the soft-clipped / not soft-clipped sequences
       revert : if true return sequences after clipping
    """
    outSeq = []
    if revert == False:
        for  (start,end) in coords:
            outSeq.append(read.seq[start:end])
    else:
        s = read.seq
        for  (start,end) in coords:
            s = s[:start] + "Z" + s[end:]
        outSeq=re.compile('Z').split(s.strip('Z'))            
    return outSeq


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("-l", "--minLen", help="Minimum length of trimmed sequence", default=0)
    parser.add_argument("-b", "--bed", help="Export soft-clipped borders coverage in BEDGraph format", action='store_true')
    parser.add_argument("-m", "--mqc", help="Export soft-clipped borders coverage in MultiQC format", action='store_true')
    parser.add_argument("-s", "--stranded", help="Annotate 3' and 5' trimmed sequences and report them in distinct output files", action='store_true')
    parser.add_argument("-v", "--verbose", help="Verbose mode", action='store_true')
    parser.add_argument("-o", "--outputDir", help="Output directory", default="./")
    args = parser.parse_args()

    if args.bed is True and args.mqc is True:
        print("--bed and --mqc parameters cannot be used together. Please specify only one coverage format")
        sys.exit(1)
    
    # Verbose mode
    if args.verbose:
        print("## extractSoftClipped.py")
        print("## BAM input=", args.filename)
        print("## minLen=", args.minLen)
        print("## outputDir=", args.outputDir)

    baseReadsFile = os.path.basename(args.filename)
    baseReadsFile = re.sub(r'\.bam$|\.sam$', '', baseReadsFile)
    ReadsFileforgenotype = re.sub(r'\.bam$|\.sam$|\.R', '', baseReadsFile)
    print('## base reads file=', baseReadsFile)
    genotype = re.search(r'^\w+-(.*)$', ReadsFileforgenotype, re.M)
    print('## genotype=', genotype.group(1))
    cur_geno = genotype.group(1)

    # Open handlers for output files
    faHandler = open(args.outputDir + '/' + baseReadsFile + '.fa', 'w')
    bkpHandler = open(args.outputDir + '/' + baseReadsFile + '_bkp.csv', 'w')
    bkpHandler.write("Read_Name" + "\t" + "Sequence" + "\t" + "Pos_Ref" + "\t" + "CigarString" + "\t" + "Label" + "\t" + "SoftClip_Position" + "\t" + "BreakPoint_Position" + "\t" +  "Breakpoint" + "\t" + "Coordinate_Fasta" +"\t" + "Length_Seq" + "\t" + "Flanking_Seq" + "\n")

    if args.bed:
        if args.stranded:
            outHandler3 = open(args.outputDir + '/' + baseReadsFile + '_3prime_bkp.bed', 'w')
            outHandler5 = open(args.outputDir + '/' + baseReadsFile + '_5prime_bkp.bed', 'w')
        else:
            outHandler = open(args.outputDir + '/' + baseReadsFile + '_bkp.bed', 'w')
    elif args.mqc:
        if args.stranded:
            outHandler3 = open(args.outputDir + '/' + baseReadsFile + '_3prime_bkp.mqc', 'w')
            outHandler5 = open(args.outputDir + '/' + baseReadsFile + '_5prime_bkp.mqc', 'w')
        else:
            outHandler = open(args.outputDir + '/' + baseReadsFile + '_bkp.mqc', 'w')


     # Read the SAM/BAM file
    if args.verbose:
        print("## Opening SAM/BAM file '", args.filename, "'...")
    samfile = pysam.Samfile(args.filename, "rb")

    # Reads are 0-based too (for both SAM and BAM format)
    # Loop on all reads
    readsCounter=0
    bkpPos = {}
 
    for read in samfile.fetch(until_eof=True):
        readsCounter =+ 1

        ## /!\ Note that here, we do not take care of read pairs
        if isSoftClipped(read):
            (clippedCoord, clippedAnnot) = getSoftClippedCoord(read, annot=args.stranded)
            clippedSeq = getClippedSeq(read, clippedCoord)

            ## Export in csv format 
            for i in range(len(clippedCoord)):
                s = clippedCoord[i]
                label = clippedAnnot[i]
                cseq = clippedSeq[i]

                if len(cseq) > int(args.minLen):
                    if label == '3prime' and args.stranded:
                        bkpStart = clippedCoord[i][0]+1
                    elif label == '5prime' and args.stranded:
                        bkpStart = clippedCoord[i][1]+1
                    else:
                        if clippedCoord[i][0]==0:
                            bkpStart = clippedCoord[i][1]+1
                        else:
                            bkpStart = clippedCoord[i][0]+1

                    ## Get bkp motif
                    mStart =  bkpStart - (int(args.minLen) + 1)
                    mEnd = bkpStart + (int(args.minLen) - 1)
                    #mStart = bkpStart - 21
                    #mEnd = bkpStart + 19
                    if mStart < 0:
                        mStart = 0
                    if mEnd >= len(read.seq):
                        mEnd = len(read.seq)
                    pos=[(mStart,mEnd)] ## 40 pb
                    bkpMotif = getClippedSeq(read, pos)
                    #if len(bkpMotif) > 0 and pos[0][1] < len(read.seq):
                    #    if len(bkpMotif[0]) > 0 :
                    #clippedSeq = getClippedSeq(read, [s])
                    if label == '5prime' and args.stranded:
                        RevCom = ReverseComplement(bkpMotif[0])
                        RevComSeq = ReverseComplement(cseq)
                        bkpHandler.write(read.query_name+ "\t" + RevCom + "\t" + str(read.reference_start+1) 
                                         + "\t" +read.cigarstring + "\t" + label + "\t" + str(s) + "\t" 
                                         + str(pos[0]) + "\t" + str(s[1])+ "\t" + str(s) + "\t" + str(s[1]) 
                                         + "\t" + RevComSeq + " \n")
                    elif (label == '3prime' and args.stranded) or not args.stranded:
                        bkpHandler.write(read.query_name+ "\t" + bkpMotif[0] +  "\t" + str(read.reference_end) 
                                         + "\t" + read.cigarstring + "\t" + label + "\t" + str(s) + "\t" 
                                         + str(pos[0]) + "\t" + str(s[0]) + "\t" + str(((s[0]),(len(read.seq)))) + "\t" 
                                         + str(len(read.seq)-s[0]) + "\t" + cseq + "\n")

            ## Export in fasta format
            ##clippedSeq = getClippedSeq(read, clippedCoord) 
            for i in range(len(clippedSeq)):
                s = clippedSeq[i]
                label = clippedAnnot[i]
                if len(s) > int(args.minLen):
                    if label == '3prime' and args.stranded:
                        faHandler.write(">" + read.query_name + "|" + str(cur_geno) + "|right|" + str(read.reference_end) + "\n")
                        faHandler.write(s + "\n")
                    elif label == '5prime' and args.stranded:
                        faHandler.write(">" + read.query_name + "|" + str(cur_geno) + "|left|" + str(read.reference_start+1) +"\n")
                        faHandler.write(s + "\n")
                    else:
                        faHandler.write(">" + read.query_name + "\n")
                        faHandler.write(s + "\n")


            coordRef = getGenomicBkp(read, excludeBorders=True, reflen=samfile.get_reference_length(samfile.get_reference_name(read.tid)))
            if args.bed or args.mqc:
                if coordRef is not None:
                    for i in range(len(coordRef)):
                        coord = coordRef[i]
                        label = clippedAnnot[i]
                        s = clippedSeq[i]

                        if len(s) > int(args.minLen):
                            if label not in bkpPos:
                                bkpPos[label] = {}
                        
                            if read.reference_name in bkpPos[label] :
                                if coord in  bkpPos[label][read.reference_name]:
                                    bkpPos[label][read.reference_name][coord] += 1
                                else:
                                    bkpPos[label][read.reference_name][coord] = 1 
                            else:
                                bkpPos[label][read.reference_name] = {}
                                bkpPos[label][read.reference_name][coord] = 1
    
        if (readsCounter % 100000 == 0 and args.verbose):
           print("##", readsCounter)
   
    ## Export BED file of breakpoints position
    if args.bed or args.mqc:
        for lab in bkpPos:
            for chr in bkpPos[lab]:
                for pos in bkpPos[lab][chr]:
                    if args.bed:
                        if args.stranded and lab == '3prime':
                            outHandler3.write(chr + "\t" + str(pos) + "\t" + str(pos+1) + "\t" + str(bkpPos[lab][chr][pos]) + "\n")
                        elif args.stranded and lab == '5prime':
                            outHandler5.write(chr + "\t" + str(pos) + "\t" + str(pos+1) + "\t" + str(bkpPos[lab][chr][pos]) + "\n")
                        else:
                            outHandler.write(chr + "\t" + str(pos) + "\t" + str(pos+1) + "\t" + str(bkpPos[lab][chr][pos]) + "\n")
                    elif args.mqc:
                        if args.stranded and lab == '3prime':
                            outHandler3.write(str(pos) + "\t" + str(bkpPos[lab][chr][pos]) + "\n")
                        elif args.stranded and lab == '5prime':
                            outHandler5.write(str(pos) + "\t" + str(bkpPos[lab][chr][pos]) + "\n")
                        else:
                            outHandler.write(str(pos) + "\t" + str(bkpPos[chr][pos]) + "\n")
    
    # Close handler
    bkpHandler.close()
    samfile.close()
    faHandler.close()
 
    if args.bed or args.mqc:
        if args.stranded:
            outHandler3.close()
            outHandler5.close()
        else:
            outHandler.close()

