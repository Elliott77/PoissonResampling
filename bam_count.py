#!/usr/bin/env python
## bam_count.py

## count reads for each allele based on known distinguishing SNPs.

import pysam
import csv
import os
import re
import sys
import time
import math
useage="""
bam_count.py <BAM> <SNP file> <out file>
"""


snp_quality_minimum = 0
##if len(sys.argv) < 3:
##    print useage
snpFileName = sys.argv[2]
original_bam = sys.argv[1]
#snp_quality_minimum = int(sys.argv[2])
out_file = sys.argv[3]

##methods
def remove_comma(snp):
    p = re.compile(r"\,")
    m = p.search(snp)
    if m:
        snp2 = snp[:m.start()]        
    else:
        snp2 = snp
    return snp2    
    
def isC57(seq, snpC57, index):
    snp_len = len(snpC57)
    if seq[index: index + snp_len] == snpC57.upper():
        return True
    else:
        return False
    
def isCAST(seq, snpCAST, index):
    snp_len = len(snpCAST)
    if seq[index: index + snp_len] == snpCAST.upper():
        return True
    else:
        return False

##snpFileName = "/uufs/chpc.utah.edu/common/home/u0368716/ref_seq/FinalVariantReferenceTable_FilterHomVarArcLiverDR_NearIndel.txt"
    
## open files
try:
    snpLines = csv.reader(open(snpFileName, "r"), delimiter ="\t")
except:
    print "error opening the file", snpFileName
    
try:
    samfile = pysam.Samfile(original_bam, "rb")
except:
    print "error opening the bam file", original_bam

f_out = open(out_file, "w")

##print "files opened"
## start loop
start_time = time.clock()    
##print "processing ", original_bam

last_chr = 'chr0'
lastSNPdepth = 5
snp_count = 0
total_read_count = 0
column_name = os.path.basename(original_bam)[0:-4]

## Column Names
#print "ChrPos\tMuscle_%s_C57\tMuscle_%s_Cast" %(column_name, column_name)
print >> f_out, "ChrPos\tMuscle_%s_C57\tMuscle_%s_Cast" %(column_name, column_name)

snpLines.next()
low_q_read = 0
low_q_snps = 0
low_q_read_indel = 0
already_written = []
aw_count = 400
last_pos = 0

for line in snpLines:

    #print"processing a SNP"
    if float(line[8]) < snp_quality_minimum:  ## skip low quality SNPs

        continue
    if line[9] == "True":  ## skip any SNPs or Indels Near Indels
        #print "skipping Varient Near Indel or Indel"
        continue
    countC57 = 0
    countCAST = 0
    count_neither = 0
    snp_count +=1
    chromosome = line[0]
    pos = long(line[1])
    chr_pos = "%s_%d" %(chromosome, pos)
    c57_snp = remove_comma(line[3])
    cast_snp = remove_comma(line[4])
    is_indel = "FALSE"
    if len(c57_snp) + len(cast_snp) > 2:
        is_indel = "TRUE"
    if chromosome != last_chr:
        elapsed_time = time.clock() - start_time


    last_chr = line[0]
    total_depth = 0
    lastSNPdepth = aw_count
    ## cull already_written list
    if len(already_written) > (2 * lastSNPdepth + 20):
        #print "Culling already written list"
        already_written = already_written[-(2 * lastSNPdepth + 1):]
    aw_count = 0
    
    for aligned_read in samfile.fetch(chromosome, int(pos)-1, int(pos)):
        aw_count += 1
        #print "checking with already written list"
        if abs(last_pos - pos) < 60 and aligned_read.qname in already_written:   ## toggle to trun on/off already written filter
            already_written.append(aligned_read.qname)
            #print "already written"
            continue  ## i.e don't count the same read twice

        total_read_count += 1                

        pos_list = aligned_read.positions
        if long(pos) - 1 in pos_list:
            index = aligned_read.positions.index(long(pos)-1)
            
            ix1 = False; l = False
            qstring = aligned_read.query   
            for i in aligned_read.cigar:    ##extract cigar data               
                if i[0]== 0:
                    ix1= i[1]
                if i[0]== 1:
                    l = i[1]
                if ix1 and l:
                    if index > ix1:    ##only if insertion is upstream of SNP
                        qstring = aligned_read.query[:ix1] + aligned_read.query[ix1+l:]   ##excize insertion
                        break
                    else:
                        qstring = aligned_read.query     ##keep the sting as is
                        break        
        else:
            continue  ## pos -1 not in list
        
        total_depth += 1

        
        if len(c57_snp) == len(cast_snp):
            is_indel = False

            if isC57(qstring, c57_snp, index):
                countC57 += 1
                already_written.append(aligned_read.qname)  
            elif isCAST(qstring, cast_snp, index):
                countCAST += 1
                already_written.append(aligned_read.qname)                            
            else:
                count_neither += 1
                already_written.append(aligned_read.qname)
                
        elif len(c57_snp) > len(cast_snp):
            aend = aligned_read.aend + 1
            pos_plus = long(pos) + len(c57_snp)            
            if aend <= pos_plus:    ## is the 3' end of the read < 3' end of the insertion
                continue            
            if isC57(qstring, c57_snp, index):
                countC57 += 1
                already_written.append(aligned_read.qname)
            elif isCAST(qstring, cast_snp, index):
                countCAST += 1
                already_written.append(aligned_read.qname)
            else:
                count_neither += 1
                already_written.append(aligned_read.qname)

        elif len(c57_snp) < len(cast_snp):
            aend = aligned_read.aend + 1
            pos_plus = long(pos) + len(cast_snp)           
            if aend <= pos_plus:    ## is the 3' end of the read < 3' end of the insertion
                count_neither += 1
                continue
            if isCAST(qstring, cast_snp, index):
                countCAST += 1
                already_written.append(aligned_read.qname)

            elif isC57(qstring, c57_snp, index):
                countC57 += 1
                already_written.append(aligned_read.qname)
            else:
                count_neither += 1
                already_written.append(aligned_read.qname)

    last_pos = pos            
    print >> f_out, "%s\t%d\t%d" %(line[6], countC57, countCAST)            
                
elapsed_time = time.clock() - start_time
print "elapsed time: %f hours\n" %(elapsed_time/3600)
samfile.close()

f_out.close()
